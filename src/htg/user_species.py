"""The user species database: the user's own species, layered over the shipped data/.

Folder (default ``platformdirs.user_data_dir("htg-equilibrium")/species``; override
with ``--user-data DIR`` or $HTG_USER_DATA, where DIR holds the files directly):

    user_pcsaft.csv     PC-SAFT rows, same columns as data/PROJECT_PCSAFT.csv
    user_epcsaft.csv    ePC-SAFT rows, same columns as data/PROJECT_EPCSAFT.csv
    user_thermo.txt     NASA-7 records (CHEMKIN 4-line format)
    user_supcrt.json    aqueous HKF and mineral Maier-Kelley records (same fields as supcrt92.xml)
    user_species.toml   [species.<name>] index (kind, formula, source), [nasa_names],
                        [solid.<name>], [[binary_interaction]], [removed]

``Model.load(user_data=DIR)`` reads data/ first, then these files
(``apply_user_layer``). A user record with the same species name replaces the
shipped one; the shipped files are never written.

The GUI uses ``UserDatabase`` (save_species, remove_species, reset_species),
``SpeciesSpec`` (the content of the new/edit species dialog), ``spec_from_model``
(fills the dialog for an existing species) and ``check_species`` (the live check list).
"""
import copy
import datetime
import json
import math
import os
import re
import tomllib
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .databases import read_nasa_lines
from .formulas import is_plain_formula, parse_formula, split_charge
from .saft import (first_complete_hkf_set, first_complete_ion_epcsaft_set, normalize_ion_name,
                   read_hkf_candidate)

APP_NAME = 'htg-equilibrium'
USER_FILES = {'pcsaft': 'user_pcsaft.csv', 'epcsaft': 'user_epcsaft.csv',
              'nasa': 'user_thermo.txt', 'supcrt': 'user_supcrt.json',
              'species': 'user_species.toml'}
RUN_RECORD_FILE = 'user_species_used.toml'

PCSAFT_COLUMNS = [
    'record_type', 'Species', 'Formula', 'Type', 'Charge', 'MW_g_mol', 'Model', 'm', 'sigma_A',
    'u_over_kB_K', 'Association_scheme', 'epsilonAB_over_kB_K', 'kappaAB', 'Dipole_D',
    'Dipole_fraction', 'Quadrupole_DA', 'Quadrupole_fraction', 'Ionization_energy_eV',
    'Source_extra_value', 'Species2', 'kij', 'Temperature dependence', 'T_unit', 'K1', 'K2',
    'RefID', 'Authors', 'Title', 'Journal', 'Year', 'DOI', 'Notes', 'Source_table']
EPCSAFT_COLUMNS = [
    'record_type', 'Species', 'Formula', 'Type', 'Charge', 'Species2', 'Model', 'Parameter',
    'Value', 'Units', 'kij', 'Temperature dependence', 'lij', 'kij_hb',
    'Dielectric correlation/value', 'Variable/type', 'RefID', 'Authors', 'Title', 'Journal',
    'Year', 'DOI', 'Notes']

# SUPCRT92 fields and units (as in supcrt92.xml; calories).
HKF_FIELDS = ['Gf', 'Hf', 'Sr', 'a1', 'a2', 'a3', 'a4', 'c1', 'c2', 'wref']
MAIER_KELLEY_FIELDS = ['Gf', 'Hf', 'Sr', 'Vr', 'a', 'b', 'c', 'Tmax']
SUPCRT_UNITS = {'Gf': 'cal/mol', 'Hf': 'cal/mol', 'Sr': 'cal/(mol*K)', 'a1': 'cal/(mol*bar)',
                'a2': 'cal/mol', 'a3': '(cal*K)/(mol*bar)', 'a4': '(cal*K)/mol',
                'c1': 'cal/(mol*K)', 'c2': '(cal*K)/mol', 'wref': 'cal/mol', 'Vr': 'cm3/mol',
                'Tmax': 'K', 'a': 'cal/(mol*K)', 'b': 'cal/(mol*K^2)', 'c': '(cal*K)/mol'}
FIELD_LABELS = {'Gf': 'ΔGf°', 'Hf': 'ΔHf°', 'Sr': 'S°', 'wref': 'ω', 'Vr': 'V°',
                'Tmax': 'maximum temperature', 'm': 'm', 'sigma': 'σ', 'epsilon': 'ε/k',
                'epsilonAB': 'εAB/k', 'kappaAB': 'κAB'}

KINDS = ('molecule', 'ion', 'solid')
KIND_LABELS = {'molecule': 'Molecule', 'ion': 'Aqueous ion', 'solid': 'Solid'}
CAN_BE_IN = {'molecule': ['gas', 'aqueous', 'organic_liquid'], 'ion': ['aqueous'], 'solid': ['solid']}
PHASE_LABELS = {'gas': 'gas', 'aqueous': 'aqueous', 'organic_liquid': 'organic liquid',
                'solid': 'solid'}
THERMO_MODES = {
    'molecule': ['nasa_record', 'nasa'],
    'ion': ['hkf_record', 'hkf'],
    'solid': ['nasa', 'maier_kelley', 'nasa_record'],
}
ASSOCIATION_SCHEMES = ['-', '2B', '3B', '4C']


# ---------------------------------------------------------------- folder
def default_user_dir():
    """$HTG_USER_DATA if set, else platformdirs.user_data_dir('htg-equilibrium')/species."""
    env = os.environ.get('HTG_USER_DATA')
    if env:
        return Path(env).expanduser()
    import platformdirs
    return Path(platformdirs.user_data_dir(APP_NAME)) / 'species'


def resolve_user_dir(path=None, enabled=True):
    """The user database folder for a command: None when switched off (--no-user-species)."""
    if not enabled:
        return None
    return Path(path).expanduser() if path else default_user_dir()


# ---------------------------------------------------------------- NASA-7 records
def _finite(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def format_nasa_record(name, elements, phase, T_low, T_mid, T_high, high, low):
    """One NASA-7 record in CHEMKIN format (4 lines of 80 columns, ending with a newline).

    high, low  the 7 coefficients of the upper (T_mid..T_high) and lower range.
    At most four elements fit in the header.
    """
    if len(elements) > 4:
        raise ValueError('A NASA-7 record holds at most 4 elements.')
    header_elements = ''
    for symbol, count in elements.items():
        count = float(count)
        text = f'{int(count):d}' if count.is_integer() else f'{count:g}'
        if len(text) > 3:
            raise ValueError(f'Element count {count:g} does not fit in a NASA-7 record.')
        header_elements += f'{symbol.upper():<2}{text:>3}'
    line1 = (f'{str(name)[:18]:<18}{"USER":<6}{header_elements:<20}{(phase or "G")[:1]:1}'
             f'{float(T_low):10.3f}{float(T_high):10.3f}{float(T_mid):8.2f}{"":6}1')
    c = [float(v) for v in high] + [float(v) for v in low]
    if len(c) != 14:
        raise ValueError('A NASA-7 record needs 7 high-range and 7 low-range coefficients.')
    fmt = '{:15.8E}'.format
    line2 = ''.join(fmt(v) for v in c[0:5]) + '    2'
    line3 = ''.join(fmt(v) for v in c[5:10]) + '    3'
    line4 = ''.join(fmt(v) for v in c[10:14]) + ' ' * 15 + '    4'
    for line in (line1, line2, line3, line4):
        if len(line) != 80:
            raise ValueError('A NASA-7 value is too large for the CHEMKIN columns.')
    return '\n'.join([line1, line2, line3, line4]) + '\n'


_NUMBER = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?')


def parse_nasa_record(text):
    """Parse one pasted NASA-7 record (4 lines, CHEMKIN format).

    Returns {name, elements, phase, T_low, T_mid, T_high, high, low}. Accepts
    lines without the 1-4 markers in column 80 and coefficients written without
    spaces between them. Raises ValueError with a plain message.
    """
    lines = [line.rstrip('\r\n') for line in str(text).splitlines()]
    lines = [line for line in lines if line.strip() and not line.strip().upper().startswith(
        ('THERMO', 'END'))]
    strict = read_nasa_lines([line + '\n' for line in lines])
    if len(strict) == 1:
        name, record = next(iter(strict.items()))
        return {'name': name, 'elements': dict(record['elements']), 'phase': record['phase'],
                'T_low': record['T_low'], 'T_mid': record['T_mid'], 'T_high': record['T_high'],
                'high': [float(v) for v in record['high']], 'low': [float(v) for v in record['low']]}
    if len(strict) > 1:
        raise ValueError('Paste one NASA-7 record (4 lines), not several.')
    if len(lines) != 4:
        raise ValueError('A NASA-7 record has 4 lines.')
    from .databases import parse_nasa_header_elements
    header = lines[0].ljust(80)
    name = header[:18].strip()
    if not name:
        raise ValueError('The first line of the NASA-7 record has no species name.')
    try:
        T_low = float(header[45:55])
        T_high = float(header[55:65])
        T_mid = float(header[65:73]) if header[65:73].strip() else 1000.0
    except ValueError:
        raise ValueError('Could not read the temperatures in the first line '
                         '(columns 46-73: T low, T high, T common).') from None
    values = []
    for line, expected in zip(lines[1:], (5, 5, 4)):
        found = [float(v.replace('D', 'E').replace('d', 'e'))
                 for v in _NUMBER.findall(line.ljust(80)[:75])]
        if len(found) < expected:
            raise ValueError(f'Line {lines.index(line) + 1} of the NASA-7 record needs '
                             f'{expected} coefficients; found {len(found)}.')
        values += found[:expected]
    return {'name': name, 'elements': parse_nasa_header_elements(header),
            'phase': header[44:45].strip().upper(), 'T_low': T_low, 'T_mid': T_mid,
            'T_high': T_high, 'high': values[0:7], 'low': values[7:14]}


def nasa_record_from_db(nasa_db, record_name):
    record = nasa_db[record_name]
    return {'name': record_name, 'elements': dict(record['elements']), 'phase': record['phase'],
            'T_low': record['T_low'], 'T_mid': record['T_mid'], 'T_high': record['T_high'],
            'high': [float(v) for v in record['high']], 'low': [float(v) for v in record['low']]}


def _split_nasa_text(text):
    """{record name: record text} of a CHEMKIN thermo file (strict columns)."""
    lines = text.splitlines(keepends=True)
    out, i = {}, 0
    while i + 3 < len(lines):
        group = [line.rstrip('\r\n') for line in lines[i:i + 4]]
        if all(len(line) >= 80 and line[79] == str(k + 1) for k, line in enumerate(group)):
            out[group[0][:18].strip()] = '\n'.join(group) + '\n'
            i += 4
        else:
            i += 1
    return out


def nasa_missing(nasa):
    """Labels of missing or invalid NASA-7 values in a coefficient dict; [] when complete."""
    nasa = nasa or {}
    missing = [label for key, label in (('T_low', 'T low'), ('T_mid', 'T common'),
                                         ('T_high', 'T high')) if not _finite(nasa.get(key))]
    for key, label in (('high', 'high-range'), ('low', 'low-range')):
        values = list(nasa.get(key) or [])
        for i in range(7):
            if i >= len(values) or not _finite(values[i]):
                missing.append(f'{label} a{i + 1}')
    return missing


# ---------------------------------------------------------------- SUPCRT records
def _elements_text(elements):
    return ''.join(f'{e}({float(n):g})' for e, n in elements.items())


def _elements_from_text(text):
    return {e: float(v) for e, v in re.findall(r'([A-Z][a-z]?)\(([-+\d.eE]+)\)', text or '')}


def supcrt_node(record):
    """XML <Species> element (as in supcrt92.xml) from a user_supcrt.json record."""
    node = ET.Element('Species')

    def sub(parent, tag, text, units=None):
        child = ET.SubElement(parent, tag)
        child.text = str(text)
        if units:
            child.set('units', units)
        return child
    sub(node, 'Name', record['Name'])
    sub(node, 'Formula', record.get('Formula') or record['Name'])
    elements = record.get('Elements') or {}
    sub(node, 'Elements', elements if isinstance(elements, str) else _elements_text(elements))
    kind = record.get('Type', 'Aqueous')
    sub(node, 'Type', kind)
    if kind == 'Aqueous':
        sub(node, 'Charge', int(record.get('Charge') or 0))
    hkf = ET.SubElement(ET.SubElement(node, 'Thermo'), 'HKF')
    sub(hkf, 'References', record.get('References', 'user'))
    sub(hkf, 'Date', record.get('Date', ''))
    values = record.get('HKF') or {}
    fields = MAIER_KELLEY_FIELDS if kind == 'Mineral' else HKF_FIELDS
    for name in fields:
        value = values.get(name)
        if value is None or not _finite(value):
            continue
        sub(hkf, name, float(value), SUPCRT_UNITS[name])
    if kind == 'Mineral':
        sub(hkf, 'NumPhaseTrans', int(values.get('NumPhaseTrans') or 0))
    return node


def supcrt_record_from_node(node):
    """user_supcrt.json record from an XML <Species> element."""
    kind = node.findtext('Type') or ''
    hkf = node.find('./Thermo/HKF')
    values = {}
    if hkf is not None:
        for child in hkf:
            if child.tag in ('References', 'Date'):
                continue
            try:
                values[child.tag] = float(child.text)
            except (TypeError, ValueError):
                values[child.tag] = None
    record = {'Name': node.findtext('Name'), 'Formula': node.findtext('Formula'), 'Type': kind,
              'Elements': _elements_from_text(node.findtext('Elements'))}
    if kind == 'Aqueous':
        record['Charge'] = int(float(node.findtext('Charge') or 0))
    record['HKF'] = values
    record['References'] = (hkf.findtext('References') if hkf is not None else '') or ''
    record['Date'] = (hkf.findtext('Date') if hkf is not None else '') or ''
    return record


def find_supcrt_node(root, name):
    return next((node for node in root.iter('Species') if node.findtext('Name') == name), None)


# ---------------------------------------------------------------- the database
def _empty_frame(columns):
    return pd.DataFrame({c: pd.Series(dtype=object) for c in columns})


def _key(name):
    return normalize_ion_name(name)


def _today():
    return datetime.date.today().isoformat()


class UserDatabase:
    """The user's species files in one folder. Nothing is written until ``save``.

    The edit methods (save_species, remove_species, ...) save immediately.
    """

    def __init__(self, directory):
        self.directory = Path(directory)
        self.pcsaft = _empty_frame(PCSAFT_COLUMNS)
        self.epcsaft = _empty_frame(EPCSAFT_COLUMNS)
        self.nasa = {}      # record name -> record text (4 lines)
        self.supcrt = []    # user_supcrt.json records
        self.config = {}    # user_species.toml
        self._normalize_config()

    def _normalize_config(self):
        c = self.config
        c.setdefault('species', {})
        c.setdefault('nasa_names', {})
        c.setdefault('solid', {})
        c.setdefault('binary_interaction', [])
        c.setdefault('removed', {}).setdefault('species', [])

    # ------------------------------------------------------------ files
    def path(self, key):
        return self.directory / USER_FILES[key]

    @classmethod
    def load(cls, directory):
        """Read the folder; missing files count as empty (the folder need not exist)."""
        if isinstance(directory, UserDatabase):
            return directory
        db = cls(directory)
        if db.path('pcsaft').is_file():
            db.pcsaft = _read_csv(db.path('pcsaft'), PCSAFT_COLUMNS)
        if db.path('epcsaft').is_file():
            db.epcsaft = _read_csv(db.path('epcsaft'), EPCSAFT_COLUMNS)
        if db.path('nasa').is_file():
            db.nasa = _split_nasa_text(db.path('nasa').read_text(encoding='utf-8', errors='ignore'))
        if db.path('supcrt').is_file():
            text = db.path('supcrt').read_text(encoding='utf-8')
            data = json.loads(text) if text.strip() else {}
            db.supcrt = list(data.get('species', []))
        if db.path('species').is_file():
            with open(db.path('species'), 'rb') as f:
                db.config = tomllib.load(f)
            db._normalize_config()
        return db

    def save(self):
        self.directory.mkdir(parents=True, exist_ok=True)
        self.pcsaft.reindex(columns=PCSAFT_COLUMNS).to_csv(self.path('pcsaft'), index=False)
        self.epcsaft.reindex(columns=EPCSAFT_COLUMNS).to_csv(self.path('epcsaft'), index=False)
        text = ''.join(self.nasa[name] for name in self.nasa)
        self.path('nasa').write_text('THERMO\n   300.000  1000.000  5000.000\n' + text + 'END\n',
                                     encoding='utf-8')
        units = {k: SUPCRT_UNITS[k] for k in dict.fromkeys(HKF_FIELDS + MAIER_KELLEY_FIELDS)}
        self.path('supcrt').write_text(json.dumps(
            {'format': 'HTG user SUPCRT92 records, version 1',
             'note': 'Fields and units as in supcrt92.xml. Type Aqueous: HKF parameters; '
                     'Type Mineral: Maier-Kelley parameters (single temperature range).',
             'units': units, 'species': self.supcrt}, indent=2, ensure_ascii=False) + '\n',
            encoding='utf-8')
        self.path('species').write_text(_dump_config(self.config), encoding='utf-8')

    # ------------------------------------------------------------ what is in it
    def owners(self):
        """{species name: set of user files with records for it}."""
        out = {}

        def add(name, key):
            if name is None or (isinstance(name, float) and math.isnan(name)) or str(name) == '':
                return
            out.setdefault(str(name), set()).add(USER_FILES[key])
        record_owner = {entry.get('nasa_record'): name
                        for name, entry in self.config['species'].items() if entry.get('nasa_record')}
        for _, row in self.pcsaft.iterrows():
            if str(row.get('record_type')) == 'pure_component':
                add(_molecule_name(row.Species), 'pcsaft')
        for _, row in self.epcsaft.iterrows():
            if str(row.get('record_type')) == 'species_parameter':
                charged = _number(row.get('Charge')) not in (None, 0.0)
                add(str(row.Species).strip() if charged else _molecule_name(row.Species), 'epcsaft')
        for record in self.nasa:
            add(record_owner.get(record, record), 'nasa')
        for record in self.supcrt:
            add(record.get('Name'), 'supcrt')
        for name in [*self.config['species'], *self.config['nasa_names'], *self.config['solid']]:
            add(name, 'species')
        return out

    def species_names(self):
        return sorted(self.owners(), key=str.casefold)

    def hidden(self):
        return list(self.config['removed']['species'])

    def has(self, name):
        return name in self.owners()

    # ------------------------------------------------------------ editing
    def remove_species(self, name, save=True):
        """Delete every user record of ``name`` (for a shipped species: back to the original)."""
        key = _key(name)
        entry = self.config['species'].pop(name, {}) or {}
        pc = self.pcsaft
        self.pcsaft = pc[~(pc.record_type.astype(str).eq('pure_component')
                           & pc.Species.map(_key).eq(key))].reset_index(drop=True)
        epc = self.epcsaft
        self.epcsaft = epc[~(epc.record_type.astype(str).eq('species_parameter')
                             & epc.Species.map(_key).eq(key))].reset_index(drop=True)
        for record in {entry.get('nasa_record'), str(name), str(name)[:18]}:
            if record and record in self.nasa and not self._nasa_record_used_by_other(record, name):
                self.nasa.pop(record)
        self.supcrt = [r for r in self.supcrt if r.get('Name') != name]
        self.config['nasa_names'].pop(name, None)
        self.config['solid'].pop(name, None)
        self.config['binary_interaction'] = [
            item for item in self.config['binary_interaction']
            if name not in item.get('species', [])]
        if save:
            self.save()

    reset_species = remove_species  # "Reset to original" for an overridden shipped species

    def _nasa_record_used_by_other(self, record, name):
        for other, entry in self.config['species'].items():
            if other != name and entry.get('nasa_record') == record:
                return True
        return any(other != name and target == record
                   for other, target in self.config['nasa_names'].items())

    def hide_species(self, name, hidden=True):
        """Hide a shipped species (listed in [removed]); hidden=False shows it again."""
        names = [n for n in self.config['removed']['species'] if n != name]
        if hidden:
            names.append(name)
        self.config['removed']['species'] = names
        self.save()

    def save_species(self, spec, model=None):
        """Store ``spec`` (a SpeciesSpec); replaces earlier user records of the same name.

        Incomplete species are saved too; they show as not usable with the reason.
        ``model`` resolves references to existing NASA-7 or SUPCRT92 records and
        the formulas of kij partners. Returns the stored species name.
        """
        spec = normalize_spec(spec)
        name = spec.name
        if not name:
            raise ValueError('Enter a name.')
        if spec.kind not in KINDS:
            raise ValueError(f'Unknown kind {spec.kind!r}.')
        try:
            elements, _ = parse_formula(spec.formula, spec.charge if spec.kind == 'ion' else None)
        except ValueError:
            elements = None
        entry = {'kind': spec.kind, 'formula': spec.formula, 'thermo': spec.thermo,
                 'updated': _today()}
        if spec.kind == 'ion':
            entry['charge'] = int(spec.charge)
        if spec.kind != 'solid':
            entry['model'] = spec.model
        for key in ('authors', 'year', 'doi', 'notes'):
            if getattr(spec, key):
                entry[key] = str(getattr(spec, key))
        keep_nasa = None
        nasa_text = None
        nasa_db = model.databases.nasa if model is not None else {}

        # Thermochemistry.
        if spec.thermo == 'nasa_record' and spec.nasa_record:
            if spec.nasa_record in self.nasa:
                keep_nasa = (spec.nasa_record, self.nasa[spec.nasa_record])
            entry['nasa_record_ref'] = spec.nasa_record
        elif spec.thermo == 'nasa':
            record_name = self._nasa_record_name(spec, nasa_db)
            entry['nasa_record'] = record_name
            missing = nasa_missing(spec.nasa)
            if missing or not elements or len(elements) > 4:
                entry['nasa_draft'] = _clean_nasa(spec.nasa)
            else:
                n = spec.nasa
                nasa_text = format_nasa_record(record_name, elements,
                                               'S' if spec.kind == 'solid' else 'G',
                                               n['T_low'], n['T_mid'], n['T_high'], n['high'], n['low'])
        supcrt_record = None
        if spec.thermo in ('hkf', 'hkf_record', 'maier_kelley'):
            values = dict(spec.maier_kelley if spec.thermo == 'maier_kelley' else spec.hkf)
            references = 'user'
            if spec.thermo == 'hkf_record' and spec.supcrt_record:
                node = find_supcrt_node(model.databases.supcrt, spec.supcrt_record) \
                    if model is not None else None
                if node is not None:
                    copied = supcrt_record_from_node(node)
                    values = {k: copied['HKF'].get(k) for k in HKF_FIELDS}
                    references = copied.get('References') or references
                entry['supcrt_record_ref'] = spec.supcrt_record
            fields = MAIER_KELLEY_FIELDS if spec.thermo == 'maier_kelley' else HKF_FIELDS
            supcrt_record = {
                'Name': name, 'Formula': spec.formula or name,
                'Type': 'Mineral' if spec.thermo == 'maier_kelley' else 'Aqueous',
                'Elements': dict(elements or {}),
                'HKF': {k: (float(values[k]) if _finite(values.get(k)) else None) for k in fields},
                'References': references, 'Date': _today()}
            if supcrt_record['Type'] == 'Aqueous':
                supcrt_record['Charge'] = int(spec.charge)
            if spec.thermo == 'hkf_record' and spec.supcrt_record:
                supcrt_record['based_on'] = spec.supcrt_record

        # Replace earlier records.
        self.remove_species(name, save=False)
        if keep_nasa:
            self.nasa[keep_nasa[0]] = keep_nasa[1]
        if nasa_text:
            self.nasa[entry['nasa_record']] = nasa_text
        if supcrt_record:
            self.supcrt.append(supcrt_record)
        record_for_name = entry.get('nasa_record') or entry.get('nasa_record_ref')
        if spec.kind == 'molecule' and record_for_name:
            self.config['nasa_names'][name] = record_for_name
        if spec.kind == 'solid' and spec.thermo in ('nasa', 'nasa_record') and record_for_name:
            self.config['solid'][name] = {'formula': spec.formula or name, 'nasa': record_for_name,
                                          'elements': dict(elements or {})}

        # Equation of state.
        if spec.kind == 'molecule' and spec.model == 'PC-SAFT':
            self.pcsaft = _append(self.pcsaft, [_pcsaft_row(spec, elements)], PCSAFT_COLUMNS)
        elif spec.kind in ('molecule', 'ion'):
            self.epcsaft = _append(self.epcsaft, _epcsaft_rows(spec), EPCSAFT_COLUMNS)
        if spec.kind == 'molecule':
            for other, kij in spec.kij:
                if not other or not _finite(kij):
                    continue
                item = {'components': [spec.formula, _eos_formula(model, other)],
                        'kij': float(kij), 'species': [name, other]}
                self.config['binary_interaction'].append(item)
        self.config['species'][name] = entry
        self.save()
        return name

    def _nasa_record_name(self, spec, nasa_db):
        """Record name for entered coefficients: the record being edited, else a free name."""
        if spec.nasa_record:
            return spec.nasa_record[:18]
        old = self.config['species'].get(spec.name, {}).get('nasa_record')
        if old:
            return old
        base = re.sub(r'\s+', '_', spec.name)[:18]
        taken = set(nasa_db) | set(self.nasa)
        candidate, n = base, 1
        while candidate in taken:
            suffix = f'_{n}'
            candidate, n = base[:18 - len(suffix)] + suffix, n + 1
        return candidate

    # ------------------------------------------------------------ low-level adds (import)
    def add_pcsaft_rows(self, rows):
        """Add pure-component rows (dicts with the shipped columns), replacing same names."""
        for row in rows:
            key = _key(row.get('Species'))
            pc = self.pcsaft
            self.pcsaft = pc[~(pc.record_type.astype(str).eq('pure_component')
                               & pc.Species.map(_key).eq(key))]
        self.pcsaft = _append(self.pcsaft, rows, PCSAFT_COLUMNS)

    def add_epcsaft_rows(self, rows):
        for key in {_key(row.get('Species')) for row in rows}:
            epc = self.epcsaft
            self.epcsaft = epc[~(epc.record_type.astype(str).eq('species_parameter')
                                 & epc.Species.map(_key).eq(key))]
        self.epcsaft = _append(self.epcsaft, rows, EPCSAFT_COLUMNS)

    def add_nasa_text(self, name, text):
        self.nasa[name] = text

    def add_supcrt_record(self, record):
        self.supcrt = [r for r in self.supcrt if r.get('Name') != record['Name']] + [record]

    def note_species(self, name, **entry):
        current = self.config['species'].get(name, {})
        self.config['species'][name] = {**current, **{k: v for k, v in entry.items() if v not in
                                                      (None, '')}, 'updated': _today()}

    # ------------------------------------------------------------ provenance
    def records_for(self, names, elements_of=None):
        """Full user records of the species ``names`` (for user_species_used.toml).

        elements_of  optional {name: elements}; user NASA-7 records with the same
                     elemental formula are included (the record a molecule uses can
                     have another name).
        """
        owners = self.owners()
        out = {}
        nasa_parsed = {name: rec for text in self.nasa.values()
                       for name, rec in read_nasa_lines(text.splitlines(keepends=True)).items()}
        wanted_formulas = [dict(e) for e in (elements_of or {}).values() if e]
        for name in dict.fromkeys(names):
            if name not in owners:
                continue
            record = {'files': sorted(owners[name])}
            entry = self.config['species'].get(name)
            if entry:
                record['index'] = copy.deepcopy(entry)
            key = _key(name)
            pc = self.pcsaft[self.pcsaft.Species.map(_key).eq(key)]
            if len(pc):
                record['pcsaft'] = _plain_rows(pc)
            epc = self.epcsaft[self.epcsaft.Species.map(_key).eq(key)]
            if len(epc):
                record['epcsaft'] = _plain_rows(epc)
            for target in [name, (entry or {}).get('nasa_record'), (entry or {}).get('nasa_record_ref'),
                           self.config['nasa_names'].get(name),
                           self.config['solid'].get(name, {}).get('nasa')]:
                if target and target in self.nasa:
                    record.setdefault('nasa', {})[target] = _nasa_plain(target, self.nasa[target],
                                                                         nasa_parsed)
            supcrt = [r for r in self.supcrt if r.get('Name') == name]
            if supcrt:
                record['supcrt'] = _drop_none(supcrt[0])
            if name in self.config['nasa_names']:
                record['nasa_name'] = self.config['nasa_names'][name]
            if name in self.config['solid']:
                record['solid'] = copy.deepcopy(self.config['solid'][name])
            kij = [item for item in self.config['binary_interaction'] if name in item.get('species', [])]
            if kij:
                record['binary_interaction'] = copy.deepcopy(kij)
            out[name] = record
        for record_name, parsed in nasa_parsed.items():
            elements = {k: float(v) for k, v in parsed['elements'].items()}
            if any({k: float(v) for k, v in e.items()} == elements for e in wanted_formulas):
                if not any(record_name in r.get('nasa', {}) for r in out.values()):
                    out.setdefault(record_name, {'files': [USER_FILES['nasa']]}).setdefault(
                        'nasa', {})[record_name] = _nasa_plain(record_name, self.nasa[record_name],
                                                               nasa_parsed)
        return out


def _molecule_name(value):
    return str(value).strip().casefold()


def _number(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _read_csv(path, columns):
    frame = pd.read_csv(path, dtype=object, keep_default_na=False, na_values=[''])
    frame.columns = [str(c).lstrip('﻿') for c in frame.columns]
    return frame.reindex(columns=list(dict.fromkeys([*columns, *frame.columns])))


def _append(frame, rows, columns):
    if not rows:
        return frame.reset_index(drop=True)
    new = pd.DataFrame(rows, dtype=object).reindex(columns=columns)
    parts = [part.astype(object) for part in (frame.reindex(columns=columns), new) if len(part)]
    return pd.concat(parts, ignore_index=True) if parts else _empty_frame(columns)


def _plain(value):
    if value is None:
        return None
    if isinstance(value, (np.floating, float)):
        return None if not math.isfinite(float(value)) else float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_plain(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    return value


def _drop_none(value):
    if isinstance(value, dict):
        return {k: _drop_none(v) for k, v in value.items() if v is not None}
    if isinstance(value, list):
        return [_drop_none(v) for v in value if v is not None]
    return value


def _plain_rows(frame):
    rows = []
    for row in frame.to_dict('records'):
        rows.append({k: _plain(v) for k, v in row.items() if _plain(v) not in (None, '')})
    return rows


def _nasa_plain(name, text, parsed):
    rec = parsed.get(name, {})
    out = {'text': text}
    if rec:
        out.update({'T_low': rec['T_low'], 'T_mid': rec['T_mid'], 'T_high': rec['T_high'],
                    'phase': rec['phase'], 'elements': {k: float(v) for k, v in rec['elements'].items()},
                    'high': [float(v) for v in rec['high']], 'low': [float(v) for v in rec['low']]})
    return out


def _clean_nasa(nasa):
    nasa = nasa or {}
    out = {}
    for key in ('T_low', 'T_mid', 'T_high'):
        if _finite(nasa.get(key)):
            out[key] = float(nasa[key])
    for key in ('high', 'low'):
        out[key] = [float(v) if _finite(v) else 'missing' for v in list(nasa.get(key) or [])][:7]
    return out


def _dump_config(config):
    import tomlkit
    doc = tomlkit.document()
    doc.add(tomlkit.comment('User species database (see docs/adding_species.md). '
                            'Written by htg.'))
    data = _drop_none(_plain(copy.deepcopy(config)))
    for key in ('species', 'nasa_names', 'solid', 'removed'):
        doc[key] = data.get(key, {})
    items = data.get('binary_interaction', [])
    if items:
        array = tomlkit.aot()
        for item in items:
            array.append(tomlkit.item(item))
        doc['binary_interaction'] = array
    for key, value in data.items():
        if key not in doc:
            doc[key] = value
    return tomlkit.dumps(doc)


def _eos_formula(model, other):
    if model is None:
        return other
    registry = model.registry
    return registry.saft_formula.get(other) or registry.saft_formula.get(str(other).casefold()) or other


def _pcsaft_row(spec, elements):
    from .formulas import molar_mass
    return {
        'record_type': 'pure_component', 'Species': spec.name, 'Formula': spec.formula,
        'Type': 'user', 'Charge': 0,
        'MW_g_mol': round(molar_mass(elements), 4) if elements else None,
        'Model': 'PC-SAFT', 'm': spec.m, 'sigma_A': spec.sigma, 'u_over_kB_K': spec.epsilon,
        'Association_scheme': spec.association_scheme or '-',
        'epsilonAB_over_kB_K': spec.epsilonAB, 'kappaAB': spec.kappaAB,
        'RefID': 'USER', 'Authors': spec.authors or None, 'Year': spec.year or None,
        'DOI': spec.doi or None, 'Notes': spec.notes or None, 'Source_table': 'user database'}


def _epcsaft_rows(spec):
    base = {'record_type': 'species_parameter', 'Species': spec.name, 'Formula': spec.formula,
            'Type': ('cation' if spec.charge > 0 else 'anion') if spec.kind == 'ion' else 'user',
            'Charge': int(spec.charge) if spec.kind == 'ion' else 0, 'Model': 'ePC-SAFT (user)',
            'RefID': 'USER', 'Authors': spec.authors or None, 'Year': spec.year or None,
            'DOI': spec.doi or None, 'Notes': spec.notes or None}
    values = [('mseg', spec.m, '-'), ('sigma', spec.sigma, 'Å'), ('u/kB', spec.epsilon, 'K')]
    if spec.kind == 'molecule':
        if spec.association_scheme and spec.association_scheme != '-':
            values.append(('association_scheme', spec.association_scheme, '-'))
        values += [('epsilonAB/kB', spec.epsilonAB, 'K'), ('kappaAB', spec.kappaAB, '-')]
    elif _finite(spec.dielectric_constant):
        values.append(('dielectric_constant', spec.dielectric_constant, '-'))
    rows = []
    for parameter, value, units in values:
        if value is None or (parameter != 'association_scheme' and not _finite(value)):
            continue
        rows.append({**base, 'Parameter': parameter, 'Value': value, 'Units': units})
    return rows


# ---------------------------------------------------------------- merging into the model
def apply_user_layer(databases, config, user_db):
    """Merge the user records into ``databases`` and the species.toml ``config`` (in place).

    Returns {species name: {'files': [...], 'overrides': bool}}: what the user
    database contributes, and whether it replaces shipped records.
    """
    owners = user_db.owners()
    shipped_pc = databases.pcsaft
    shipped_epc = databases.epcsaft
    shipped_nasa = set(databases.nasa)
    shipped_supcrt = {node.findtext('Name') for node in databases.supcrt.iter('Species')}
    shipped_pc_names = set(shipped_pc[shipped_pc.record_type.astype(str).eq('pure_component')]
                           .Species.map(_key))
    shipped_epc_names = set(shipped_epc[shipped_epc.record_type.astype(str).eq('species_parameter')]
                            .Species.map(_key))
    shipped_solids = set(config.get('solid', {}))

    # SAFT tables: user rows replace shipped rows of the same species and come first,
    # so they win where the lookup goes by formula.
    pc_user = user_db.pcsaft[user_db.pcsaft.record_type.astype(str).eq('pure_component')]
    epc_user = user_db.epcsaft[user_db.epcsaft.record_type.astype(str).eq('species_parameter')]
    names = set(pc_user.Species.map(_key)) | set(epc_user.Species.map(_key))
    if len(pc_user):
        databases.pcsaft = _merge_frame(shipped_pc, pc_user, 'pure_component', names,
                                        USER_FILES['pcsaft'])
    if len(epc_user) or len(pc_user):
        databases.epcsaft = _merge_frame(shipped_epc, epc_user, 'species_parameter', names,
                                         USER_FILES['epcsaft'])

    # NASA-7: same record name replaces.
    for text in user_db.nasa.values():
        databases.nasa.update(read_nasa_lines(text.splitlines(keepends=True)))

    # SUPCRT92: user records replace same names and go first.
    if user_db.supcrt:
        root = databases.supcrt
        user_names = {r.get('Name') for r in user_db.supcrt}
        for node in [n for n in root.findall('Species') if n.findtext('Name') in user_names]:
            root.remove(node)
        children = list(root)
        position = next((i for i, n in enumerate(children) if n.tag == 'Species'), len(children))
        for offset, record in enumerate(user_db.supcrt):
            root.insert(position + offset, supcrt_node(record))

    # species.toml tables.
    nasa_names = config.setdefault('nasa_names', {})
    for name, record in user_db.config['nasa_names'].items():
        if record in databases.nasa:
            nasa_names[name] = record
    solids = config.setdefault('solid', {})
    for name, entry in user_db.config['solid'].items():
        if {'formula', 'nasa', 'elements'} <= set(entry):
            solids[name] = dict(entry)
    pairs = config.setdefault('binary_interaction', [])
    for item in user_db.config['binary_interaction']:
        components = item.get('components') or []
        if len(components) != 2 or not _finite(item.get('kij')):
            continue
        pair = set(components)
        pairs[:] = [p for p in pairs if set(p['components']) != pair]
        pairs.append({'components': list(components), 'kij': float(item['kij'])})

    info = {}
    for name, files in owners.items():
        key = _key(name)
        overrides = (key in shipped_pc_names or key in shipped_epc_names or name in shipped_nasa
                     or name in shipped_supcrt or name in shipped_solids)
        entry = user_db.config['species'].get(name, {})
        for record in (entry.get('nasa_record'),):
            overrides = overrides or (record in shipped_nasa)
        info[name] = {'files': sorted(files), 'overrides': bool(overrides)}
    return info


def _merge_frame(shipped, user, record_type, names, source_file):
    drop = shipped.record_type.astype(str).eq(record_type) & shipped.Species.map(_key).isin(names)
    kept = shipped[~drop]
    user = user.copy()
    user['source_file'] = source_file
    user.index = range(-len(user), 0)
    columns = list(dict.fromkeys([*shipped.columns, 'source_file']))
    user = user.reindex(columns=columns)
    for column in shipped.columns:
        if pd.api.types.is_numeric_dtype(shipped[column]):
            user[column] = pd.to_numeric(user[column], errors='coerce')
    merged = pd.concat([user.astype(object), kept.reindex(columns=columns).astype(object)])
    for column in shipped.columns:
        if pd.api.types.is_numeric_dtype(shipped[column]):
            merged[column] = pd.to_numeric(merged[column], errors='coerce')
    return merged


def hide_species(registry, names):
    """Drop hidden shipped species from a discovered registry and its catalogue."""
    names = set(names)
    if not names:
        return
    for table in (registry.molecular, registry.organic, registry.solid, registry.ion):
        for name in names & set(table):
            table.pop(name)
    registry.catalogue = registry.catalogue[~registry.catalogue.Species.isin(names)] \
        .reset_index(drop=True)


# ---------------------------------------------------------------- the dialog content
@dataclass
class SpeciesSpec:
    """Content of the new/edit species dialog.

    kind     'molecule', 'ion' or 'solid'
    thermo   molecule: 'nasa_record' (use nasa_record) or 'nasa' (enter coefficients)
             ion:      'hkf_record' (copy supcrt_record) or 'hkf' (enter hkf)
             solid:    'nasa' (phase S coefficients), 'maier_kelley' or 'nasa_record'
    nasa     {'T_low', 'T_mid', 'T_high' [K], 'high': 7 coefficients, 'low': 7 coefficients}
    hkf      {Gf, Hf, Sr, a1, a2, a3, a4, c1, c2, wref} in SUPCRT92 units (cal, bar)
    maier_kelley  {Gf, Hf, Sr, Vr, a, b, c, Tmax} in SUPCRT92 units (cal, cm3, K)
    model    'PC-SAFT' or 'ePC-SAFT' (ions: ePC-SAFT only)
    m, sigma [Å], epsilon [K], association_scheme, epsilonAB [K], kappaAB
    kij      [(other species name, kij)]
    """
    kind: str = 'molecule'
    name: str = ''
    formula: str = ''
    charge: int = 0
    thermo: str = 'nasa_record'
    nasa_record: str = ''
    nasa: dict = field(default_factory=dict)
    supcrt_record: str = ''
    hkf: dict = field(default_factory=dict)
    maier_kelley: dict = field(default_factory=dict)
    model: str = 'PC-SAFT'
    m: float | None = None
    sigma: float | None = None
    epsilon: float | None = None
    association_scheme: str = '-'
    epsilonAB: float | None = None
    kappaAB: float | None = None
    dielectric_constant: float | None = None
    kij: list = field(default_factory=list)
    authors: str = ''
    year: str = ''
    doi: str = ''
    notes: str = ''

    def to_dict(self):
        return asdict(self)


def normalize_spec(spec):
    """A copy with trimmed text; molecule names are lower case (as the registry stores them)."""
    if isinstance(spec, dict):
        spec = SpeciesSpec(**spec)
    spec = copy.deepcopy(spec)
    spec.name = str(spec.name or '').strip()
    if spec.kind == 'molecule':
        spec.name = spec.name.casefold()
    spec.formula = str(spec.formula or '').strip()
    spec.charge = int(spec.charge or 0)
    if spec.kind == 'ion':
        spec.model = 'ePC-SAFT'
    spec.year = str(spec.year or '').strip()
    for key in ('m', 'sigma', 'epsilon', 'epsilonAB', 'kappaAB', 'dielectric_constant'):
        value = getattr(spec, key)
        setattr(spec, key, float(value) if _finite(value) else None)
    spec.association_scheme = str(spec.association_scheme or '-').strip() or '-'
    return spec


def spec_from_model(model, name, reference_temperature_C=None):
    """SpeciesSpec filled from what ``model`` uses for ``name`` (to edit an existing species)."""
    registry, databases = model.registry, model.databases
    T = float(reference_temperature_C if reference_temperature_C is not None
              else model.load_args.get('reference_temperature_C') or 350.0) + 273.15
    user_entry = {}
    if model.user_data is not None:
        user_entry = UserDatabase.load(model.user_data).config['species'].get(name, {})
    spec = SpeciesSpec(name=name)
    if name in registry.ion or name in registry.ion_charge:
        spec.kind, spec.model = 'ion', 'ePC-SAFT'
        spec.charge = int(registry.ion_charge.get(name, registry.ion.get(name, {}).get('charge', 0)))
        spec.formula = registry.ion.get(name, {}).get('formula') or registry.ion_saft_formula.get(name, name)
        epc = first_complete_ion_epcsaft_set(databases.epcsaft, name, registry.ion_saft_formula.get(name),
                                             spec.charge)
        if epc:
            spec.m, spec.sigma, spec.epsilon = epc['m'], epc['sigma'], epc['epsilon']
            spec.dielectric_constant = epc['dielectric_constant'] if _finite(
                epc['dielectric_constant']) else None
            _source_from_row(spec, databases.epcsaft, epc['Database row'])
        hkf = first_complete_hkf_set(databases.supcrt, name, registry.supcrt_ion_aliases.get(name, []))
        if hkf:
            spec.thermo, spec.supcrt_record = 'hkf_record', hkf['SUPCRT name']
            spec.hkf = {k: hkf[k] for k in HKF_FIELDS}
        else:
            spec.thermo = 'hkf'
    elif name in registry.saft_formula or (name in registry.all_molecular and name not in registry.solid):
        spec.kind = 'molecule'
        entry = registry.all_molecular.get(name, {})
        spec.formula = registry.saft_formula.get(name) or entry.get('formula') or name
        params = registry.saft_parameters(databases, name, T) if name in registry.saft_formula else None
        if params:
            spec.model = 'ePC-SAFT' if 'EPCSAFT' in str(params['Source']).upper() else 'PC-SAFT'
            spec.m, spec.sigma, spec.epsilon = params['m'], params['sigma'], params['epsilon']
            spec.association_scheme = params.get('association_scheme') or '-'
            spec.epsilonAB = params['epsilonAB'] if _finite(params['epsilonAB']) else None
            spec.kappaAB = params['kappaAB'] if _finite(params['kappaAB']) else None
            frame = databases.epcsaft if spec.model == 'ePC-SAFT' else databases.pcsaft
            _source_from_row(spec, frame, params['Database row'])
            if name == 'water' and spec.model == 'ePC-SAFT':
                spec.notes = ('σ is temperature dependent in the shipped data; shown at the '
                              'reference temperature.')
        record = entry.get('nasa')
        if record is None:
            rows = registry.catalogue[(registry.catalogue.Species == name)
                                      & (registry.catalogue.Phase == 'molecular')]
            record = rows.iloc[0].NASA if len(rows) else None
        if record is None or record not in databases.nasa:
            record = {str(k).casefold(): v for k, v in registry.config.get('nasa_names', {}).items()} \
                .get(name, record)
        if record in databases.nasa:
            spec.thermo, spec.nasa_record = 'nasa_record', record
            spec.nasa = _nasa_dict(nasa_record_from_db(databases.nasa, record))
        else:
            spec.thermo = 'nasa'
        formulas = {}
        for other in [*registry.all_molecular, *registry.saft_formula]:
            formulas.setdefault(registry.saft_formula.get(other), other)
        for (a, b), kij in registry.binary_interactions:
            if spec.formula in (a, b):
                other = b if a == spec.formula else a
                spec.kij.append((formulas.get(other, other), kij))
    elif name in registry.solid:
        spec.kind = 'solid'
        entry = registry.solid[name]
        spec.formula = entry.get('formula') or name
        if entry.get('adapter') == 'supcrt_mineral':
            spec.thermo = 'maier_kelley'
            spec.maier_kelley = dict(entry['supcrt_mineral'])
        else:
            record = entry.get('nasa')
            spec.thermo, spec.nasa_record = 'nasa_record', record or ''
            if record in databases.nasa:
                spec.nasa = _nasa_dict(nasa_record_from_db(databases.nasa, record))
    else:
        row = registry.catalogue[registry.catalogue.Species == name]
        if not len(row):
            raise KeyError(f'Unknown species {name!r}.')
        spec.formula = str(row.iloc[0].Formula)
        phase = str(row.iloc[0].Phase)
        spec.kind = 'solid' if 'Mineral' in phase or phase == 'solid' else (
            'ion' if 'Aqueous' in phase or phase == 'aqueous ion' else 'molecule')
        if spec.kind == 'molecule':
            spec.thermo = 'nasa_record' if row.iloc[0].get('NASA') in databases.nasa else 'nasa'
            spec.nasa_record = row.iloc[0].get('NASA') or '' if spec.thermo == 'nasa_record' else ''
            params = registry.saft_parameters(databases, name, T) if name in registry.saft_formula else None
            if params:
                spec.m, spec.sigma, spec.epsilon = params['m'], params['sigma'], params['epsilon']
        elif spec.kind == 'ion':
            node = find_supcrt_node(databases.supcrt, name)
            if node is not None:
                spec.charge = int(float(node.findtext('Charge') or 0))
                spec.thermo, spec.supcrt_record = 'hkf_record', name
            spec.model = 'ePC-SAFT'
    if user_entry:
        for key in ('authors', 'year', 'doi', 'notes'):
            if user_entry.get(key):
                setattr(spec, key, str(user_entry[key]))
        if user_entry.get('nasa_draft') and spec.thermo != 'nasa_record':
            spec.thermo, spec.nasa = 'nasa', dict(user_entry['nasa_draft'])
    return spec


def _nasa_dict(record):
    return {k: record[k] for k in ('T_low', 'T_mid', 'T_high', 'high', 'low')}


def _source_from_row(spec, frame, index):
    try:
        row = frame.loc[index]
    except KeyError:
        return
    for key, column in (('authors', 'Authors'), ('year', 'Year'), ('doi', 'DOI')):
        value = row.get(column)
        if value is not None and not (isinstance(value, float) and math.isnan(value)) and str(value) != 'nan':
            text = str(value)
            if key == 'year' and text.endswith('.0'):
                text = text[:-2]
            setattr(spec, key, text)


# ---------------------------------------------------------------- live check
@dataclass
class CheckItem:
    ok: bool
    text: str


@dataclass
class CheckResult:
    items: list
    can_be_in: list
    usable: bool

    @property
    def can_be_in_text(self):
        return 'Can be in: ' + ', '.join(PHASE_LABELS[p] for p in self.can_be_in)


def nasa_records_for_formula(model, formula, solid=False):
    """NASA-7 record names with the elemental formula of ``formula`` (in database order).

    solid  False: gas/liquid records (for molecules); True: phase S records; None: all.
    """
    try:
        elements, _ = parse_formula(formula)
    except ValueError:
        return []
    target = {k: float(v) for k, v in elements.items()}
    return [name for name, rec in model.databases.nasa.items()
            if {k: float(v) for k, v in rec['elements'].items()} == target
            and (solid is None or (rec['phase'] == 'S') == solid)]


def nasa_match_message(model, formula, record=None, solid=False):
    """'Only record with formula CH4O', '3 records with formula C3H6O', or 'No record with formula …'.

    solid  as in nasa_records_for_formula (the dialog passes kind == 'solid').
    """
    from .formulas import hill_formula
    try:
        elements, _ = parse_formula(formula)
    except ValueError as exc:
        return str(exc)
    text = hill_formula(elements)
    names = nasa_records_for_formula(model, formula, solid)
    if record and record not in names:
        if record in model.databases.nasa:
            other = hill_formula(model.databases.nasa[record]['elements'])
            return f'Record {record} has formula {other}, not {text}'
        return f'No NASA-7 record named {record}'
    if not names:
        return f'No record with formula {text}'
    if len(names) == 1:
        return f'Only record with formula {text}'
    return f'{len(names)} records with formula {text}'


def supcrt_aqueous_records(model, charge=None):
    """Names of complete aqueous HKF records in supcrt92.xml (with the given charge)."""
    out = []
    for node in model.databases.supcrt.iter('Species'):
        if node.findtext('Type') != 'Aqueous':
            continue
        values = read_hkf_candidate(node)
        if values is None:
            continue
        if charge is not None and int(values['Charge']) != int(charge):
            continue
        out.append(node.findtext('Name'))
    return out


def _missing_text(prefix, labels, limit=6):
    shown = ', '.join(labels[:limit]) + (f' and {len(labels) - limit} more' if len(labels) > limit else '')
    return f'{prefix}: {shown}'


def check_species(spec, model, editing=None):
    """The live check list of the species dialog.

    editing  name of the species being edited (shipped or user); the name check
             then accepts that name.
    """
    spec = normalize_spec(spec)
    items = []
    # 1. Formula.
    elements = None
    try:
        elements, _ = parse_formula(spec.formula, spec.charge if spec.kind == 'ion' else None)
        if spec.kind == 'ion' and not spec.charge:
            items.append(CheckItem(False, 'Enter the charge of the ion.'))
        elif spec.kind == 'molecule' and not is_plain_formula(spec.formula):
            items.append(CheckItem(False, 'Write the formula without brackets or dots, for example '
                                          'C2H6O. The equation of state uses it as the component name.'))
        else:
            items.append(CheckItem(True, 'Formula understood'))
    except ValueError as exc:
        implied = split_charge(spec.formula)[1] if spec.kind == 'ion' else None
        if implied is not None and spec.charge and implied != spec.charge:
            items.append(CheckItem(False, f'The charge in the formula ({implied:+d}) differs from '
                                          f'Charge ({spec.charge:+d}).'))
        else:
            items.append(CheckItem(False, str(exc)))

    # 2. Thermochemistry.
    items.append(_check_thermo(spec, model, elements))

    # 3. Equation of state.
    if spec.kind != 'solid':
        items.append(_check_eos(spec, model))

    # 4. Name.
    items.append(_check_name(spec, model, editing))
    return CheckResult(items=items, can_be_in=list(CAN_BE_IN.get(spec.kind, [])),
                       usable=all(item.ok for item in items))


def _check_thermo(spec, model, elements):
    nasa_db = model.databases.nasa
    ok = CheckItem(True, 'Thermochemistry complete')
    mode = spec.thermo
    if mode not in THERMO_MODES.get(spec.kind, []):
        return CheckItem(False, 'Choose how to give the thermochemistry.')
    if mode == 'nasa_record':
        if not spec.nasa_record:
            return CheckItem(False, 'Choose a NASA-7 record.')
        if spec.nasa_record not in nasa_db:
            return CheckItem(False, f'NASA-7 record {spec.nasa_record} not found.')
        record = nasa_db[spec.nasa_record]
        if spec.kind == 'solid' and record['phase'] != 'S':
            return CheckItem(False, f'NASA-7 record {spec.nasa_record} is not a solid (phase S) record.')
        if spec.kind == 'molecule' and record['phase'] == 'S':
            return CheckItem(False, f'NASA-7 record {spec.nasa_record} is a solid record.')
        if elements is not None and {k: float(v) for k, v in record['elements'].items()} != \
                {k: float(v) for k, v in elements.items()}:
            from .formulas import hill_formula
            return CheckItem(False, f'NASA-7 record {spec.nasa_record} has formula '
                                    f'{hill_formula(record["elements"])}, not {hill_formula(elements)}.')
        return ok
    if mode == 'nasa':
        missing = nasa_missing(spec.nasa)
        if missing:
            return CheckItem(False, _missing_text('Missing NASA-7 values', missing))
        n = spec.nasa
        if not float(n['T_low']) < float(n['T_mid']) < float(n['T_high']):
            return CheckItem(False, 'NASA-7 temperatures must rise: T low < T common < T high.')
        if elements is not None and len(elements) > 4:
            return CheckItem(False, 'A NASA-7 record holds at most 4 elements.')
        if elements is None:
            return CheckItem(False, 'Thermochemistry needs a formula that is understood.')
        return ok
    if mode == 'hkf_record':
        if not spec.supcrt_record:
            return CheckItem(False, 'Choose a SUPCRT92 aqueous record.')
        node = find_supcrt_node(model.databases.supcrt, spec.supcrt_record)
        if node is None:
            return CheckItem(False, f'SUPCRT92 record {spec.supcrt_record} not found.')
        values = read_hkf_candidate(node) if node.findtext('Type') == 'Aqueous' else None
        if values is None:
            return CheckItem(False, f'SUPCRT92 record {spec.supcrt_record} has no complete aqueous HKF data.')
        if int(values['Charge']) != int(spec.charge):
            return CheckItem(False, f'Charge of {spec.supcrt_record} ({int(values["Charge"]):+d}) '
                                    f'differs from {int(spec.charge):+d}.')
        return ok
    fields = HKF_FIELDS if mode == 'hkf' else MAIER_KELLEY_FIELDS
    values = spec.hkf if mode == 'hkf' else spec.maier_kelley
    missing = [FIELD_LABELS.get(k, k) for k in fields if not _finite((values or {}).get(k))]
    if missing:
        return CheckItem(False, _missing_text('Missing HKF values' if mode == 'hkf'
                                              else 'Missing Maier–Kelley values', missing))
    if mode == 'maier_kelley':
        if float(values['Vr']) <= 0:
            return CheckItem(False, 'V° must be positive.')
        if float(values['Tmax']) < 298.15:
            return CheckItem(False, 'The maximum temperature must be at least 298.15 K.')
    return ok


def _check_eos(spec, model):
    if spec.kind == 'ion' and spec.model != 'ePC-SAFT':
        return CheckItem(False, 'Ions need ePC-SAFT parameters.')
    missing = [FIELD_LABELS[k] for k in ('m', 'sigma', 'epsilon')
               if getattr(spec, k) is None or getattr(spec, k) <= 0]
    if missing:
        return CheckItem(False, _missing_text('Missing equation of state parameters', missing))
    if spec.kind == 'molecule' and spec.association_scheme not in ('-', ''):
        assoc = [FIELD_LABELS[k] for k in ('epsilonAB', 'kappaAB')
                 if getattr(spec, k) is None or getattr(spec, k) <= 0]
        if assoc:
            return CheckItem(False, f'Association scheme {spec.association_scheme} needs '
                                    + ' and '.join(assoc) + '.')
    if spec.kind == 'molecule':
        known = set(model.registry.all_molecular)
        for other, kij in spec.kij:
            if not other:
                continue
            if other not in known:
                return CheckItem(False, f'kij: unknown species {other}.')
            if not _finite(kij):
                return CheckItem(False, f'kij with {other} is not a number.')
    return CheckItem(True, 'Equation of state parameters complete')


def _check_name(spec, model, editing):
    if not spec.name:
        return CheckItem(False, 'Enter a name.')
    if editing is not None and _key(spec.name) == _key(editing):
        return CheckItem(True, 'Name not used by another species')
    used = {_key(n): n for n in model.registry.catalogue.Species.astype(str)}
    if _key(spec.name) in used:
        return CheckItem(False, f'Name already used by {used[_key(spec.name)]}. '
                                'Edit that species, or choose another name.')
    return CheckItem(True, 'Name not used by another species')


# ---------------------------------------------------------------- listing
def user_species_table(model):
    """Species from the user database: Species, Kind, Origin ('yours' or 'overrides shipped'), Files."""
    rows = []
    entries = {}
    if model.user_data is not None:
        entries = UserDatabase.load(model.user_data).config['species']
    for name, info in sorted(model.user_info.items(), key=lambda kv: kv[0].casefold()):
        rows.append({'Species': name, 'Kind': KIND_LABELS.get(entries.get(name, {}).get('kind'), ''),
                     'Origin': 'overrides shipped' if info['overrides'] else 'yours',
                     'Files': ', '.join(info['files'])})
    return pd.DataFrame(rows, columns=['Species', 'Kind', 'Origin', 'Files'])


# ---------------------------------------------------------------- saved runs
def used_species_names(sweep, model_args=None):
    """Species a run can use: the case's phase lists and feed species, the species of
    the solved points and (with an aqueous phase) the aqueous reference species."""
    case = sweep.case
    names = ['water']
    phases = getattr(case, 'phases', None) or {}
    for species in phases.values():
        names += list(species or [])
    feed_species = getattr(getattr(case, 'feed', None), 'species_moles', None)
    if feed_species:
        names += list(feed_species)
    elements_of = {}
    for result in sweep.results:
        for (_, species), data in result.species_data.items():
            names.append(species)
            elements_of[species] = data.get('elements', {})
    if phases.get('aqueous'):
        species_file = None
        if model_args:
            species_file = model_args.get('species_file') or (
                Path(model_args['data_dir']) / 'species.toml' if model_args.get('data_dir') else None)
        if species_file is None:
            from .model import default_data_dir
            species_file = default_data_dir() / 'species.toml'
        try:
            with open(species_file, 'rb') as f:
                reference = tomllib.load(f).get('aqueous_reference', {})
            names += list(reference.get('neutrals', [])) + list(reference.get('ions', []))
        except OSError:
            pass
    return list(dict.fromkeys(names)), elements_of


def run_record(model_args, sweep):
    """What a saved run records about the user database; None when it was not used.

    model_args  Model.load arguments (Model.load_args), with 'user_data'.
    """
    user_data = (model_args or {}).get('user_data')
    if not user_data:
        return None
    db = UserDatabase.load(user_data)
    names, elements_of = used_species_names(sweep, model_args)
    return {'database': str(Path(user_data).resolve()), 'records': db.records_for(names, elements_of)}


def write_run_record(directory, record):
    """Write user_species_used.toml (when user records were used) and note the
    user database folder at the end of case.toml."""
    import tomlkit
    directory = Path(directory)
    if record.get('records'):
        doc = tomlkit.document()
        doc.add(tomlkit.comment('User species records (full parameters) used by this run.'))
        doc['database'] = record['database']
        doc['species'] = _drop_none(_plain(record['records']))
        (directory / RUN_RECORD_FILE).write_text(tomlkit.dumps(doc), encoding='utf-8')
    case_file = directory / 'case.toml'
    if case_file.exists():
        text = case_file.read_text(encoding='utf-8')
        note = f'\n# User species database: {record["database"]}\n'
        if record.get('records'):
            note += f'# User species used by this run: see {RUN_RECORD_FILE}\n'
        case_file.write_text(text.rstrip('\n') + '\n' + note, encoding='utf-8')
