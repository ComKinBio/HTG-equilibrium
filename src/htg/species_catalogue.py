"""The Species tab table: one row per species, plain-language status and the filters.

    table = build_catalogue(model)
    shown = filter_catalogue(table, text='meth', elements='C, H, O', exact=False,
                             source='Any', show='Usable only', can_be_in=['gas'])
    count_text(shown, table)      # '6 of 412 species'

Columns of the table:
    Species, Formula, Elements ({symbol: count}), Kind ('molecule', 'aqueous ion',
    'solid' or 'record'), Can be in (list of phase keys), Data (list of database
    labels), Yours (bool: records in the user database), Overrides (bool: they
    replace shipped records), Feed (bool: usable as a feed species), Usable (bool),
    Status ('Usable' or the reason in plain language).
"""
import re

import pandas as pd

from .formulas import hill_formula, parse_element_list, parse_formula
from .saft import first_complete_ion_epcsaft_set, normalize_ion_name
from .species import NO_SAFT_REASON

PHASE_LABELS = {'gas': 'Gas', 'aqueous': 'Aqueous', 'organic_liquid': 'Organic liquid', 'solid': 'Solid'}
DATA_SOURCES = ['NASA-7', 'PC-SAFT', 'ePC-SAFT', 'SUPCRT92']
SOURCE_OPTIONS = ['Any', *DATA_SOURCES]
SHOW_OPTIONS = ['Usable only', 'Not usable only', 'All']
USABLE = 'Usable'

# Registry catalogue phase -> phases the species can be in (the old TARGETS of the GUI).
TARGETS = {'molecular': ['gas', 'aqueous', 'organic_liquid'], 'aqueous ion': ['aqueous'],
           'solid': ['solid'], 'SUPCRT Aqueous': ['aqueous'], 'SUPCRT Gaseous': ['gas'],
           'SUPCRT Mineral': ['solid']}
KINDS = {'molecular': 'molecule', 'aqueous ion': 'aqueous ion', 'solid': 'solid'}

# Registry reasons -> plain text. None: worked out from the databases (see _reason).
READY_REASONS = {'ready', 'NASA solid', 'SUPCRT mineral: single-range Maier-Kelley'}
PLAIN_REASONS = {
    NO_SAFT_REASON: 'No PC-SAFT parameters',
    'Complete SAFT parameters missing': 'No complete PC-SAFT parameters',
    'Atomic weight missing': 'Unknown element in the formula',
    'Elemental composition unresolved': 'Formula not understood',
    'NASA mapping missing/ambiguous': None,
    'Complete matching ionic EOS/HKF/element data missing': None,
    'SUPCRT reference present; matching EOS/reference adapter required': None,
    'NASA solid': None,
    'Mineral HKF data missing': 'No Maier–Kelley data',
    'Polymorphic mineral: transition adapter not implemented': 'Has phase transitions (not supported yet)',
    'Invalid/missing mineral parameters': 'Maier–Kelley data incomplete',
    'Element composition or atomic weight missing': 'Unknown element in the formula',
    'Expected Vr in cm3/mol (do not use mislabeled MolarVolume)': 'Unexpected units in the SUPCRT92 record',
}
PLAIN_PATTERNS = [
    (re.compile(r'^Unexpected units for mineral'), 'Unexpected units in the SUPCRT92 record'),
    (re.compile(r'^could not convert string to float'), 'Maier–Kelley data incomplete'),
    (re.compile(r'float\(\) argument must be'), 'Maier–Kelley data incomplete'),
]

COLUMNS = ['Species', 'Formula', 'Elements', 'Kind', 'Can be in', 'Data', 'Yours', 'Overrides',
           'Feed', 'Usable', 'Status']


# ---------------------------------------------------------------- building
def build_catalogue(model):
    """One row per species of ``model`` (ready or not); see the module docstring."""
    registry, databases = model.registry, model.databases
    nasa_db = databases.nasa
    lookup = _Lookup(model)
    rows = {}
    for record in registry.catalogue.to_dict('records'):
        name = str(record['Species'])
        phase = str(record['Phase'])
        row = rows.get(name)
        if row is None:
            row = rows[name] = {'Species': name, 'phases': [], 'ready': [], 'reasons': [],
                                'formula': record.get('Formula'), 'nasa': None}
        row['phases'].append(phase)
        row['ready'].append(bool(record['Ready']))
        nasa = record.get('NASA')
        if isinstance(nasa, str) and nasa in nasa_db:
            row['nasa'] = row['nasa'] or nasa
        if not record['Ready']:
            row['reasons'].append(_reason(model, lookup, name, phase, str(record['Reason']), record))
    feed = {*registry.all_molecular, *registry.solid}
    user = getattr(model, 'user_info', {}) or {}
    out = []
    for name, row in rows.items():
        ready_phases = [p for p, r in zip(row['phases'], row['ready']) if r]
        usable = bool(ready_phases)
        targets = []
        for phase in (ready_phases or row['phases']):
            targets += TARGETS.get(phase, [])
        targets = [p for p in PHASE_LABELS if p in targets]
        main = ready_phases[0] if ready_phases else row['phases'][0]
        kind = KINDS.get(main, 'record')
        if main == 'SUPCRT Mineral':
            kind = 'solid'
        elements = _elements(registry, lookup, name, row, nasa_db)
        formula = str(row['formula']) if row['formula'] is not None else ''
        if row['nasa'] and formula == row['nasa'] and elements:
            formula = hill_formula(elements)
        reasons = list(dict.fromkeys(r for r in row['reasons'] if r))
        info = user.get(name)
        out.append({
            'Species': name, 'Formula': formula, 'Elements': elements, 'Kind': kind,
            'Can be in': targets, 'Data': lookup.sources(name, formula, row, kind),
            'Yours': info is not None, 'Overrides': bool(info and info['overrides']),
            'Feed': usable and name in feed, 'Usable': usable,
            'Status': USABLE if usable else ('; '.join(reasons) or 'Not usable')})
    return pd.DataFrame(out, columns=COLUMNS)


class _Lookup:
    """Name and formula sets of the databases, built once per catalogue."""

    def __init__(self, model):
        databases = model.databases
        pc = databases.pcsaft
        pc = pc[pc.record_type.astype(str).eq('pure_component')]
        self.pc_names = set(pc.Species.map(normalize_ion_name))
        self.pc_formulas = set(pc.Formula.map(normalize_ion_name))
        epc = databases.epcsaft
        epc = epc[epc.record_type.astype(str).eq('species_parameter')]
        self.epc_names = set(epc.Species.map(normalize_ion_name)) | set(epc.Formula.map(normalize_ion_name))
        self.supcrt = {}
        for node in databases.supcrt.iter('Species'):
            self.supcrt.setdefault(node.findtext('Name'), node)
        self.supcrt_aliases = {alias for aliases in model.registry.supcrt_ion_aliases.values()
                               for alias in aliases}
        self.ion_aliases = model.registry.supcrt_ion_aliases
        nasa_by_formula = {}
        for name, record in databases.nasa.items():
            key = tuple(sorted((k, float(v)) for k, v in record.get('elements', {}).items()))
            nasa_by_formula.setdefault(key, []).append(name)
        self.nasa_by_formula = nasa_by_formula
        self.nasa = databases.nasa

    def sources(self, name, formula, row, kind):
        key = normalize_ion_name(name)
        found = []
        if row['nasa'] or name in self.nasa:
            found.append('NASA-7')
        if key in self.pc_names or (kind in ('molecule',) and normalize_ion_name(formula) in self.pc_formulas):
            found.append('PC-SAFT')
        if key in self.epc_names:
            found.append('ePC-SAFT')
        if name in self.supcrt or self.ion_aliases.get(name):
            found.append('SUPCRT92')
        return found


def _elements(registry, lookup, name, row, nasa_db):
    for table in (registry.ion, registry.molecular, registry.organic, registry.solid):
        if name in table and table[name].get('elements'):
            return {k: _count(v) for k, v in table[name]['elements'].items()}
    if row['nasa']:
        return {k: _count(v) for k, v in nasa_db[row['nasa']]['elements'].items()}
    node = lookup.supcrt.get(name)
    if node is not None:
        found = re.findall(r'([A-Z][a-z]?)\(([-+\d.eE]+)\)', node.findtext('Elements') or '')
        if found:
            return {e: _count(v) for e, v in found}
    try:
        return {k: _count(v) for k, v in parse_formula(row['formula'])[0].items()}
    except (ValueError, TypeError):
        return {}


def _count(value):
    value = float(value)
    return int(value) if value.is_integer() else value


def _reason(model, lookup, name, phase, reason, record):
    """Plain-language version of a registry reason (the mapping table plus database lookups)."""
    if reason in PLAIN_REASONS and PLAIN_REASONS[reason] is not None:
        return PLAIN_REASONS[reason]
    for pattern, text in PLAIN_PATTERNS:
        if pattern.search(reason):
            return text
    if reason == 'NASA mapping missing/ambiguous':
        formula = str(record.get('Formula') or '')
        try:
            elements = parse_formula(formula)[0]
        except ValueError:
            return 'Formula not understood'
        key = tuple(sorted((k, float(v)) for k, v in elements.items()))
        candidates = lookup.nasa_by_formula.get(key, [])
        text = hill_formula(elements)
        if not candidates:
            return f'No NASA-7 record with formula {text}'
        if len(candidates) > 1:
            return f'Several NASA-7 records with formula {text}; choose one'
        return f'Another species also has formula {text}; choose its NASA-7 record'
    if reason == 'Complete matching ionic EOS/HKF/element data missing':
        registry = model.registry
        charge = registry.ion_charge.get(name, 0)
        if not first_complete_ion_epcsaft_set(model.databases.epcsaft, name,
                                              registry.ion_saft_formula.get(name), charge):
            return 'No complete ePC-SAFT ion parameters'
        if not registry.supcrt_ion_aliases.get(name):
            return f'No SUPCRT92 aqueous record with charge {charge:+d}'
        return 'SUPCRT92 record incomplete'
    if reason == 'SUPCRT reference present; matching EOS/reference adapter required':
        node = lookup.supcrt.get(name)
        kind = node.findtext('Type') if node is not None else ''
        if kind == 'Aqueous':
            charge = float((node.findtext('Charge') if node is not None else 0) or 0)
            return 'No ePC-SAFT ion parameters' if charge else 'No PC-SAFT parameters'
        if kind == 'Gaseous':
            return 'SUPCRT92 gas record only; gases use NASA-7 and PC-SAFT'
        return 'Only a SUPCRT92 record'
    if reason == 'NASA solid':
        return f'NASA-7 record {record.get("NASA")} not found'
    return reason


# ---------------------------------------------------------------- filters
def filter_catalogue(table, text='', elements='', exact=False, source='Any', show='Usable only',
                     can_be_in=()):
    """Rows of ``table`` (from build_catalogue) that pass the Species tab filters.

    text       part of the name or formula (case-insensitive)
    elements   'C, H, O': species containing at least these elements; with
               exact=True only species made of exactly these elements
    source     'Any' or one of DATA_SOURCES
    show       'Usable only', 'Not usable only' or 'All'
    can_be_in  phase keys or labels ('gas', 'Organic liquid', ...); a species
               matches if it can be in any of them; empty: no filter
    Raises ValueError for unknown element symbols or options.
    """
    mask = pd.Series(True, index=table.index)
    text = str(text or '').strip().casefold()
    if text:
        mask &= (table.Species.astype(str).str.casefold().str.contains(text, regex=False)
                 | table.Formula.astype(str).str.casefold().str.contains(text, regex=False))
    symbols = set(parse_element_list(elements))
    if symbols:
        if exact:
            mask &= table.Elements.map(lambda e: bool(e) and set(e) == symbols)
        else:
            mask &= table.Elements.map(lambda e: symbols <= set(e or {}))
    if source and source != 'Any':
        if source not in DATA_SOURCES:
            raise ValueError(f'Unknown data source {source!r}.')
        mask &= table.Data.map(lambda found: source in found)
    if show not in SHOW_OPTIONS:
        raise ValueError(f'Unknown show option {show!r}.')
    if show == 'Usable only':
        mask &= table.Usable
    elif show == 'Not usable only':
        mask &= ~table.Usable
    phases = {_phase_key(p) for p in (can_be_in or [])}
    if phases:
        mask &= table['Can be in'].map(lambda targets: bool(phases & set(targets)))
    return table[mask]


def _phase_key(value):
    text = str(value).strip()
    for key, label in PHASE_LABELS.items():
        if text in (key, label) or text.casefold() == label.casefold():
            return key
    raise ValueError(f'Unknown phase {value!r}.')


def count_text(shown, total):
    """'6 of 412 species' (shown and total: tables or numbers)."""
    n = shown if isinstance(shown, int) else len(shown)
    m = total if isinstance(total, int) else len(total)
    return f'{n} of {m} species'


def table_rows(table):
    """JSON-friendly rows for a NiceGUI table: lists as text, labels for phases."""
    rows = []
    for record in table.to_dict('records'):
        rows.append({**record,
                     'Elements': ' · '.join(f'{e} {n:g}' for e, n in record['Elements'].items()),
                     'Can be in': [PHASE_LABELS[p] for p in record['Can be in']],
                     'Data': list(record['Data']) + (['yours'] if record['Yours'] else [])})
    return rows
