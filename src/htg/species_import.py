"""Import species from a file into the user database (Species tab, "Import from file…").

    fmt = detect_format(text)                       # 'pcsaft', 'epcsaft', 'nasa', 'supcrt' or None
    rows = preview_import(text, model, fmt)         # one ImportRow per species in the file
    apply_import([r for r in rows if r.selected], model.user_data)

Each preview row says what would happen: 'new', 'replaces existing values'
(with the parameters that differ) or 'skipped' (with the reason), and whether
the species would be usable after the import. Rows that replace values start
unticked; skipped rows cannot be ticked.
"""
import io
import math
import shutil
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .databases import read_nasa_lines
from .formulas import is_plain_formula, parse_formula
from .saft import (PCSAFT_MODEL_PRIORITY, first_complete_epcsaft_set, first_complete_ion_epcsaft_set,
                   normalize_ion_name, read_hkf_candidate)
from .user_species import (EPCSAFT_COLUMNS, HKF_FIELDS, MAIER_KELLEY_FIELDS, PCSAFT_COLUMNS,
                           UserDatabase, _split_nasa_text, find_supcrt_node, supcrt_record_from_node)

FORMATS = {'pcsaft': 'PC-SAFT table', 'epcsaft': 'ePC-SAFT table', 'nasa': 'NASA-7 thermo file',
           'supcrt': 'SUPCRT92 XML'}
NEW, REPLACES, SKIPPED = 'new', 'replaces existing values', 'skipped'
NOTE = ('Species that replace existing values are unticked by default. Imported species go into '
        'your own database and can be removed again on the Species tab.')

PCSAFT_FIELDS = [('Formula', 'formula'), ('m', 'm'), ('sigma_A', 'σ'), ('u_over_kB_K', 'ε/k'),
                 ('Association_scheme', 'association scheme'), ('epsilonAB_over_kB_K', 'εAB/k'),
                 ('kappaAB', 'κAB')]
EOS_LABELS = {'Formula': 'formula', 'm': 'm', 'sigma': 'σ', 'epsilon': 'ε/k',
              'association_scheme': 'association scheme', 'epsilonAB': 'εAB/k', 'kappaAB': 'κAB',
              'dielectric_constant': 'dielectric constant'}


@dataclass
class ImportRow:
    key: str
    species: str
    formula: str
    action: str                 # NEW, REPLACES or SKIPPED
    reason: str = ''            # why it is skipped
    after: str = ''             # 'Usable' or the reason, after the import
    differs: list = field(default_factory=list)
    fmt: str = ''
    payload: dict = field(default_factory=dict, repr=False)

    @property
    def selectable(self):
        return self.action != SKIPPED

    @property
    def selected(self):
        """Default tick: new species only."""
        return self.action == NEW

    @property
    def what_happens(self):
        return f'{SKIPPED}: {self.reason}' if self.action == SKIPPED else self.action

    def as_row(self):
        """Row for the preview table."""
        return {'key': self.key, 'Species': self.species, 'Formula': self.formula,
                'What happens': self.what_happens, 'After import': self.after,
                'selectable': self.selectable, 'selected': self.selected}


# ---------------------------------------------------------------- format
def detect_format(text):
    """'pcsaft', 'epcsaft', 'nasa', 'supcrt' or None, from the file content."""
    text = str(text or '').lstrip('﻿')
    head = text.lstrip()
    if head.startswith('<') and '<Species' in text:
        return 'supcrt'
    first = head.splitlines()[0] if head else ''
    columns = {c.strip().strip('"').casefold() for c in first.split(',')}
    if {'parameter', 'value', 'species'} <= columns:
        return 'epcsaft'
    if 'species' in columns and ({'m', 'sigma_a'} & columns):
        return 'pcsaft'
    if read_nasa_lines(text.splitlines(keepends=True)) or _split_nasa_text(text):
        return 'nasa'
    return None


def read_text(path):
    return Path(path).read_text(encoding='utf-8-sig', errors='replace')


# ---------------------------------------------------------------- preview
def preview_import(text, model, fmt=None):
    """ImportRow list for the file content ``text`` (format detected unless ``fmt``)."""
    fmt = fmt or detect_format(text)
    if fmt not in FORMATS:
        raise ValueError('File format not recognised. Choose the format.')
    rows = {'pcsaft': _preview_pcsaft, 'epcsaft': _preview_epcsaft, 'nasa': _preview_nasa,
            'supcrt': _preview_supcrt}[fmt](str(text).lstrip('﻿'), model)
    for row in rows:
        row.fmt = fmt
    _after_import(rows, model)
    return rows


def _number(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _text(value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ''
    return str(value).strip()


def _same(a, b):
    na, nb = _number(a), _number(b)
    if na is not None or nb is not None:
        return na is not None and nb is not None and math.isclose(na, nb, rel_tol=1e-9, abs_tol=1e-12)
    return _text(a).casefold() in ('', '-', 'nan') and _text(b).casefold() in ('', '-', 'nan') \
        or _text(a).casefold() == _text(b).casefold()


def _read_table(text):
    frame = pd.read_csv(io.StringIO(text), dtype=object, keep_default_na=False, na_values=[''])
    frame.columns = [str(c).strip().lstrip('﻿') for c in frame.columns]
    return frame


def _current_eos(model, name, T=623.15):
    """Equation-of-state parameters the model now has for molecule ``name`` (None if none)."""
    registry = model.registry
    key = normalize_ion_name(name)
    if name in registry.saft_formula:
        params = registry.saft_parameters(model.databases, name, T)
        if params:
            return params
    pc = model.databases.pcsaft
    rows = pc[pc.record_type.astype(str).eq('pure_component') & pc.Species.map(normalize_ion_name).eq(key)]
    if len(rows):
        row = rows.iloc[0]
        cols = model.databases.pc_columns
        return {'Formula': row.get('Formula'), 'm': row.get(cols['m']), 'sigma': row.get(cols['sigma']),
                'epsilon': row.get(cols['epsilon']),
                'association_scheme': row.get(cols['association']) if cols['association'] else '-',
                'epsilonAB': row.get(cols['epsilonAB']) if cols['epsilonAB'] else None,
                'kappaAB': row.get(cols['kappaAB']) if cols['kappaAB'] else None}
    epc = model.databases.epcsaft
    rows = epc[epc.record_type.astype(str).eq('species_parameter') & epc.Species.map(normalize_ion_name).eq(key)]
    if len(rows):
        return {'Formula': rows.iloc[0].get('Formula')}
    return None


def _differences(current, new):
    return [EOS_LABELS[k] for k in EOS_LABELS if k in new and not _same(current.get(k), new.get(k))]


def _preview_pcsaft(text, model):
    frame = _read_table(text)
    missing = [c for c in ('record_type', 'Species', 'Formula', 'm', 'sigma_A', 'u_over_kB_K')
               if c not in frame.columns]
    if missing:
        raise ValueError(f'The PC-SAFT table lacks the columns: {", ".join(missing)}.')
    rows, seen = [], set()
    for index, record in frame.iterrows():
        raw_name = _text(record.get('Species'))
        name = raw_name.casefold()
        formula = _text(record.get('Formula'))
        row = ImportRow(key=f'pcsaft:{index}', species=name or f'(row {index + 2})', formula=formula,
                        action=NEW)
        kind = _text(record.get('record_type'))
        reason = ''
        if kind != 'pure_component':
            reason = ('binary interaction rows are not imported; add kij in the species dialog'
                      if kind == 'binary_interaction' else f'record type "{kind}" is not imported')
        elif not name:
            reason = 'no species name'
        elif not formula:
            reason = 'no formula'
        elif not is_plain_formula(formula):
            reason = 'formula not understood (use element symbols and counts, e.g. C2H6O)'
        else:
            try:
                parse_formula(formula)
            except ValueError as exc:
                reason = str(exc).rstrip('.')
        if not reason:
            for column, label in (('m', 'm'), ('sigma_A', 'σ'), ('u_over_kB_K', 'ε/k')):
                value = _number(record.get(column))
                if value is None or value <= 0:
                    reason = f'{label} is not a positive number'
                    break
        if not reason:
            model_name = _text(record.get('Model')).casefold() or 'pc-saft'
            if model_name not in PCSAFT_MODEL_PRIORITY:
                reason = f'model "{_text(record.get("Model"))}" is not used (PC-SAFT or GC-PPC-SAFT)'
        if not reason and name in seen:
            reason = 'listed twice in the file; the first row is used'
        if reason:
            row.action, row.reason = SKIPPED, reason
            rows.append(row)
            continue
        seen.add(name)
        new = {'Formula': formula, 'm': record.get('m'), 'sigma': record.get('sigma_A'),
               'epsilon': record.get('u_over_kB_K'),
               'association_scheme': _text(record.get('Association_scheme')) or '-',
               'epsilonAB': record.get('epsilonAB_over_kB_K'), 'kappaAB': record.get('kappaAB')}
        current = _current_eos(model, name)
        if current is not None:
            row.differs = _differences(current, new)
            if row.differs:
                row.action = REPLACES
            else:
                row.action, row.reason = SKIPPED, 'same values as the existing record'
        payload = {c: (None if pd.isna(record.get(c)) else record.get(c)) for c in PCSAFT_COLUMNS
                   if c in frame.columns}
        payload['Species'] = raw_name
        payload.setdefault('Model', 'PC-SAFT')
        if not payload.get('Model'):
            payload['Model'] = 'PC-SAFT'
        row.payload = {'pcsaft': [payload], 'kind': 'molecule'}
        rows.append(row)
    return rows


def _preview_epcsaft(text, model):
    frame = _read_table(text)
    missing = [c for c in ('record_type', 'Species', 'Parameter', 'Value') if c not in frame.columns]
    if missing:
        raise ValueError(f'The ePC-SAFT table lacks the columns: {", ".join(missing)}.')
    rows = []
    species_rows = frame[frame.record_type.astype(str).eq('species_parameter')]
    other = frame[~frame.record_type.astype(str).eq('species_parameter')]
    for kind, part in other.groupby(other.record_type.astype(str), sort=False):
        rows.append(ImportRow(key=f'epcsaft:{kind}', species=f'{len(part)} {kind} row(s)', formula='',
                              action=SKIPPED, reason=f'record type "{kind}" is not imported'))
    for number, (raw_name, group) in enumerate(species_rows.groupby('Species', sort=False)):
        raw_name = _text(raw_name)
        charge = _number(group['Charge'].dropna().iloc[0]) if 'Charge' in group and \
            group['Charge'].notna().any() else 0.0
        charge = charge or 0.0
        formula = _text(group['Formula'].dropna().iloc[0]) if 'Formula' in group and \
            group['Formula'].notna().any() else ''
        name = raw_name if charge else raw_name.casefold()
        row = ImportRow(key=f'epcsaft:{number}', species=name, formula=formula, action=NEW)
        typed = group.copy()
        typed['Charge'] = pd.to_numeric(typed.get('Charge'), errors='coerce')
        if charge:
            complete = first_complete_ion_epcsaft_set(typed, raw_name, formula, int(charge))
        else:
            complete = first_complete_epcsaft_set(typed, name, formula, [name], 623.15) \
                if formula else None
        reason = ''
        if not formula:
            reason = 'no formula'
        elif complete is None:
            reason = 'incomplete ePC-SAFT parameters (need m, σ and ε/k)'
        elif not charge and not is_plain_formula(formula):
            reason = 'formula not understood (use element symbols and counts, e.g. C2H6O)'
        if reason:
            row.action, row.reason = SKIPPED, reason
            rows.append(row)
            continue
        if charge:
            exists = name in model.registry.ion_charge
            current = first_complete_ion_epcsaft_set(model.databases.epcsaft, name,
                                                     model.registry.ion_saft_formula.get(name, formula),
                                                     int(charge)) if exists else None
        else:
            current = _current_eos(model, name)
            exists = current is not None
        if exists:
            row.differs = _differences(current or {}, {k: v for k, v in complete.items()
                                                       if k in EOS_LABELS and k != 'Formula'})
            if row.differs or current is None:
                row.action = REPLACES
            else:
                row.action, row.reason = SKIPPED, 'same values as the existing record'
        payload = [{c: (None if pd.isna(r.get(c)) else r.get(c)) for c in EPCSAFT_COLUMNS if c in group}
                   for r in group.to_dict('records')]
        row.payload = {'epcsaft': payload, 'kind': 'ion' if charge else 'molecule', 'charge': charge}
        rows.append(row)
    return rows


def _preview_nasa(text, model):
    records = _split_nasa_text(text)
    if not records:
        raise ValueError('No NASA-7 records found (4 lines each, CHEMKIN format).')
    nasa_db = model.databases.nasa
    rows = []
    for number, (name, record_text) in enumerate(records.items()):
        parsed = read_nasa_lines(record_text.splitlines(keepends=True)).get(name)
        row = ImportRow(key=f'nasa:{number}', species=name, formula='', action=NEW)
        reason = ''
        if parsed is None:
            reason = 'record not understood'
        else:
            values = [parsed['T_low'], parsed['T_mid'], parsed['T_high'], *parsed['low'], *parsed['high']]
            if not np.all(np.isfinite(values)):
                reason = 'missing coefficients'
            elif not parsed['T_low'] < parsed['T_mid'] < parsed['T_high']:
                reason = 'temperatures must rise: T low < T common < T high'
            elif not parsed['elements']:
                reason = 'no elements in the header'
        if parsed is not None and parsed['elements']:
            from .formulas import hill_formula
            row.formula = hill_formula({k: (int(v) if float(v).is_integer() else v)
                                        for k, v in parsed['elements'].items()})
        if reason:
            row.action, row.reason = SKIPPED, reason
            rows.append(row)
            continue
        if name in nasa_db:
            old = nasa_db[name]
            differs = [label for key, label in (('T_low', 'T low'), ('T_mid', 'T common'),
                                                 ('T_high', 'T high'), ('phase', 'phase'))
                       if not _same(old[key], parsed[key])]
            if not np.allclose(old['high'], parsed['high'], rtol=1e-9, atol=0):
                differs.append('high-range coefficients')
            if not np.allclose(old['low'], parsed['low'], rtol=1e-9, atol=0):
                differs.append('low-range coefficients')
            if dict(old['elements']) != dict(parsed['elements']):
                differs.append('elements')
            row.differs = differs
            if differs:
                row.action = REPLACES
            else:
                row.action, row.reason = SKIPPED, 'same values as the existing record'
        row.payload = {'nasa': (name, record_text),
                       'kind': 'solid' if parsed['phase'] == 'S' else 'molecule'}
        rows.append(row)
    return rows


def _preview_supcrt(text, model):
    try:
        root = ET.fromstring(text)
    except ET.ParseError as exc:
        raise ValueError(f'The XML file could not be read: {exc}.') from None
    nodes = [root] if root.tag == 'Species' else list(root.iter('Species'))
    if not nodes:
        raise ValueError('No <Species> records found in the XML file.')
    rows = []
    for number, node in enumerate(nodes):
        record = supcrt_record_from_node(node)
        name = record['Name'] or ''
        row = ImportRow(key=f'supcrt:{number}', species=name, formula=record.get('Formula') or '',
                        action=NEW)
        kind = record['Type']
        reason = ''
        if not name:
            reason = 'no name'
        elif kind == 'Aqueous':
            if read_hkf_candidate(node) is None:
                reason = 'incomplete HKF data'
        elif kind == 'Mineral':
            values = record['HKF']
            if int(values.get('NumPhaseTrans') or 0):
                reason = 'has phase transitions (not supported yet)'
            elif any(_number(values.get(k)) is None for k in MAIER_KELLEY_FIELDS):
                reason = 'incomplete Maier–Kelley data'
        else:
            reason = f'{kind or "unknown"} records are not used'
        if reason:
            row.action, row.reason = SKIPPED, reason
            rows.append(row)
            continue
        existing = find_supcrt_node(model.databases.supcrt, name)
        if existing is not None:
            old = supcrt_record_from_node(existing)
            fields = HKF_FIELDS if kind == 'Aqueous' else MAIER_KELLEY_FIELDS
            differs = [f for f in fields if not _same(old['HKF'].get(f), record['HKF'].get(f))]
            if kind == 'Aqueous' and old.get('Charge') != record.get('Charge'):
                differs.insert(0, 'charge')
            row.differs = differs
            if differs:
                row.action = REPLACES
            else:
                row.action, row.reason = SKIPPED, 'same values as the existing record'
        row.payload = {'supcrt': record, 'kind': 'ion' if kind == 'Aqueous' else 'solid'}
        rows.append(row)
    return rows


# ---------------------------------------------------------------- after import
def _after_import(rows, model):
    """Fill ImportRow.after by loading the model with the rows applied to a copy of the user database."""
    importable = [r for r in rows if r.action != SKIPPED]
    if not importable:
        return
    from .model import Model
    from .species_catalogue import build_catalogue
    with tempfile.TemporaryDirectory(prefix='htg-import-') as tmp:
        target = Path(tmp) / 'species'
        if model.user_data is not None and Path(model.user_data).is_dir():
            shutil.copytree(model.user_data, target)
        apply_import(importable, target)
        args = {**model.load_args, 'user_data': str(target)}
        catalogue = build_catalogue(Model.load(**args)).set_index('Species')
    for row in importable:
        status = catalogue.Status.get(row.species, 'Not listed after import')
        row.after = status + (f'; differs: {", ".join(row.differs)}' if row.differs else '')


# ---------------------------------------------------------------- apply
def apply_import(rows, user_data):
    """Write the (selected) rows to the user database in ``user_data``; returns how many."""
    db = UserDatabase.load(user_data)
    count = 0
    for row in rows:
        if row.action == SKIPPED:
            continue
        payload = row.payload
        entry = {'kind': payload.get('kind'), 'formula': row.formula, 'imported': row.fmt}
        if 'pcsaft' in payload:
            db.add_pcsaft_rows(payload['pcsaft'])
            first = payload['pcsaft'][0]
            entry.update(model='PC-SAFT', authors=_text(first.get('Authors')),
                         year=_text(first.get('Year')), doi=_text(first.get('DOI')))
        if 'epcsaft' in payload:
            db.add_epcsaft_rows(payload['epcsaft'])
            first = payload['epcsaft'][0]
            entry.update(model='ePC-SAFT', authors=_text(first.get('Authors')),
                         year=_text(first.get('Year')), doi=_text(first.get('DOI')))
            if payload.get('charge'):
                entry['charge'] = int(payload['charge'])
        if 'nasa' in payload:
            name, text = payload['nasa']
            db.add_nasa_text(name, text)
            entry.update(thermo='nasa', nasa_record=name)
        if 'supcrt' in payload:
            db.add_supcrt_record(payload['supcrt'])
            record = payload['supcrt']
            entry.update(thermo='maier_kelley' if record['Type'] == 'Mineral' else 'hkf')
            if record.get('Charge'):
                entry['charge'] = int(record['Charge'])
        db.note_species(row.species, **entry)
        count += 1
    db.save()
    return count
