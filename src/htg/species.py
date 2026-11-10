"""Species registry: which species exist, their composition and database identifiers.

The registry starts from data/species.toml and is completed by scanning the
databases (``Registry.discover``). It never invents data: a species is "ready"
only when every required record exists. ``Registry.catalogue`` lists every
candidate, ready or not, with the reason.
"""
import copy
import re
import tomllib

import numpy as np
import pandas as pd

from .constants import ATOMIC_WEIGHTS
from .saft import (first_complete_epcsaft_set, first_complete_hkf_set,
                   first_complete_ion_epcsaft_set, first_complete_pcsaft_set, read_hkf_candidate)

ELEMENT_PATTERN = r'([A-Z][a-z]?)\(([-+\d.eE]+)\)'


def formula_elements(formula):
    """Elemental composition of a plain formula like 'C2H6O'; None if it cannot be parsed."""
    text = str(formula).strip()
    if not re.fullmatch(r'(?:[A-Z][a-z]?\d*)+', text):
        return None
    result = {}
    for symbol, number in re.findall(r'([A-Z][a-z]?)(\d*)', text):
        result[symbol] = result.get(symbol, 0.0) + int(number or 1)
    return result


def molecular_weight(formula, atomic_weights=ATOMIC_WEIGHTS):
    tokens = re.findall('([A-Z][a-z]?)(\\d*)', str(formula).strip())
    if not tokens:
        raise ValueError(f'Could not parse molecular formula: {formula}')
    total = 0.0
    for element, number in tokens:
        if element not in atomic_weights:
            raise KeyError(f"Atomic weight missing for '{element}' in formula '{formula}'.")
        total += atomic_weights[element] * (int(number) if number else 1)
    return float(total)


def _ion_identity(name, charge):
    text = re.sub(r'\s+', '', str(name)).replace('^', '')
    # The explicit charge avoids treating the 2 in VO2+ as a divalent charge.
    sign = '+' if charge > 0 else '-'
    suffixes = [f'{abs(int(charge))}{sign}', f'{sign}{abs(int(charge))}', sign * abs(int(charge))]
    if text.upper() in ('NA+', 'K+', 'LI+', 'RB+', 'CS+', 'CL-', 'F-', 'BR-', 'I-'):
        return text.rstrip('+-').casefold(), int(charge)
    for suffix in suffixes:
        if text.endswith(suffix):
            return text[:-len(suffix)].casefold(), int(charge)
    return text.rstrip('+-').casefold(), int(charge)


class Registry:
    """All species known to the model.

    molecular, organic  fluid molecules {name: {formula, nasa, elements}}
    solid               pure solids {name: {formula, nasa | supcrt_mineral, elements, ...}}
    ion                 aqueous ions {name: {formula, charge, supcrt, elements}}
    """

    def __init__(self, config):
        self.config = config
        self.molecular = {}
        self.organic = {}
        self.solid = {}
        self.ion = {}
        # Identifiers used to look species up in the SAFT and SUPCRT databases.
        self.saft_formula = {}
        self.saft_aliases = {}
        self.ion_charge = {}
        self.ion_saft_formula = {}
        self.supcrt_ion_aliases = {}
        for group in ('molecular', 'organic'):
            for name, entry in config.get(group, {}).items():
                getattr(self, group)[name] = {k: entry[k] for k in ('formula', 'nasa', 'elements')}
                self.saft_formula[name] = entry['formula']
                self.saft_aliases[name] = list(entry['saft_aliases'])
        for name, entry in config.get('solid', {}).items():
            self.solid[name] = {k: entry[k] for k in ('formula', 'nasa', 'elements')}
        for name, entry in config.get('ion', {}).items():
            self.ion[name] = {k: entry[k] for k in ('formula', 'charge', 'supcrt', 'elements')}
            self.ion_charge[name] = entry['charge']
            self.ion_saft_formula[name] = entry['saft_formula']
            self.supcrt_ion_aliases[name] = list(entry['supcrt_aliases'])
        reference = config['aqueous_reference']
        self.aqueous_reference_neutrals = list(reference['neutrals'])
        self.aqueous_reference_ions = list(reference['ions'])
        self.probe_weights = dict(reference['probe_weights'])
        self.binary_interactions = [(tuple(item['components']), float(item['kij']))
                                    for item in config.get('binary_interaction', [])]
        self.catalogue = pd.DataFrame()
        self.ion_epcsaft_table = pd.DataFrame()
        self.ion_hkf_table = pd.DataFrame()

    @classmethod
    def from_file(cls, path, nasa_names=None):
        with open(path, 'rb') as f:
            config = tomllib.load(f)
        if nasa_names is not None:
            config['nasa_names'] = dict(nasa_names)
        return cls(config)

    # ---------------------------------------------------------------- lookups
    @property
    def all_molecular(self):
        return {**self.molecular, **self.organic}

    def nasa_name(self, name):
        registry = {**self.molecular, **self.organic, **self.solid}
        if name not in registry:
            raise KeyError(f'No NASA mapping for {name}')
        return registry[name]['nasa']

    def entry(self, phase, name):
        """Registry entry of ``name`` in ``phase``; raises with the available names otherwise."""
        if phase == 'solid':
            registry = self.solid
        elif phase == 'aqueous' and name in self.ion:
            registry = self.ion
        else:
            registry = self.all_molecular
        if name not in registry:
            raise ValueError(f"{phase}: '{name}' is not registered. "
                             f"Available species: {', '.join(registry)}")
        return registry[name]

    def any_entry(self, name):
        return {**self.all_molecular, **self.solid, **self.ion}[name]

    def available(self, include_unavailable=False):
        if include_unavailable:
            return self.catalogue.copy()
        return self.catalogue[self.catalogue.Ready].reset_index(drop=True)

    # ---------------------------------------------------------------- SAFT parameters
    def saft_parameters(self, databases, name, T):
        """First complete ePC-SAFT set, else PC-SAFT set, for molecule ``name`` at T [K]."""
        formula = self.saft_formula.get(name)
        return (first_complete_epcsaft_set(databases.epcsaft, name, formula,
                                           self.saft_aliases.get(name), T)
                or first_complete_pcsaft_set(databases, name, formula))

    def resolve_parameters(self, databases, T, species, atomic_weights=ATOMIC_WEIGHTS):
        """SAFT parameters of every molecule in ``species`` at T [K]; raises when incomplete."""
        resolved = {}
        all_molecular = self.all_molecular
        for name in dict.fromkeys(species):
            if name not in all_molecular:
                raise ValueError(f'{name}: not found in the molecular database registry')
            data = self.saft_parameters(databases, name, T)
            if data is None:
                raise ValueError(f'{name}: NASA data alone is insufficient for PC-SAFT; '
                                 'complete SAFT parameters/mapping required')
            data['Mw'] = sum(atomic_weights[e] * v
                             for e, v in all_molecular[name]['elements'].items())
            if name == 'water' and data.get('association_scheme') in ('-', '', None):
                data['association_scheme'] = '2B'
            resolved[name] = data
        return resolved

    # ---------------------------------------------------------------- discovery
    def discover(self, databases, T):
        """Add every species with complete database records; build the catalogue.

        T [K] is used only to evaluate temperature-dependent SAFT parameters
        while checking completeness.
        """
        AW = ATOMIC_WEIGHTS
        nasa_db = databases.nasa
        overrides = {str(k).strip().casefold(): {'nasa': v}
                     for k, v in self.config.get('nasa_names', {}).items()}
        original_molecules = {**self.molecular, **self.organic}
        nasa_by_formula = {}
        for nasa_name, data in nasa_db.items():
            key = tuple(sorted(data.get('elements', {}).items()))
            if key:
                nasa_by_formula.setdefault(key, []).append(nasa_name)
        pcsaft, epcsaft = databases.pcsaft, databases.epcsaft
        pc = pcsaft[pcsaft['record_type'].astype(str).eq('pure_component')]
        epc = epcsaft[epcsaft['record_type'].astype(str).eq('species_parameter')
                      & epcsaft['Charge'].fillna(0).eq(0)]
        records = []
        for frame in (pc, epc):
            for _, row in frame.iterrows():
                if pd.isna(row.get('Species')) or pd.isna(row.get('Formula')):
                    continue
                name = str(row.Species).strip().casefold()
                if name not in {r[0] for r in records}:
                    records.append((name, str(row.Formula).strip()))

        # Fluid molecules: SAFT record + unambiguous NASA record.
        catalogue = []
        for name, formula in records:
            elements = formula_elements(formula)
            override = overrides.get(name, {})
            mapped = override.get('nasa') or original_molecules.get(name, {}).get('nasa')
            if not mapped:
                exact = next((key for key in nasa_db if key.casefold() in (name, formula.casefold())),
                             None)
                if exact:
                    mapped = exact
                elif elements:
                    candidates = nasa_by_formula.get(tuple(sorted(elements.items())), [])
                    same_formula = [other for other, f in records if formula_elements(f) == elements]
                    if len(candidates) == 1 and len(same_formula) == 1:
                        mapped = candidates[0]
            if override.get('nasa') and mapped not in nasa_db:
                raise ValueError(f'NASA override for {name} not found: {mapped}')
            if mapped in nasa_db:
                elements = nasa_db[mapped].get('elements') or elements
            if not elements:
                catalogue.append(dict(Species=name, Formula=formula, Phase='molecular', Ready=False,
                                      Reason='Elemental composition unresolved', NASA=mapped))
                continue
            self.saft_formula[name] = formula
            self.saft_aliases[name] = [name, formula]
            data = self.saft_parameters(databases, name, T)
            ready = bool(mapped in nasa_db and data and all(e in AW for e in elements))
            reason = ('ready' if ready else 'NASA mapping missing/ambiguous' if mapped not in nasa_db
                      else 'Complete SAFT parameters missing' if not data else 'Atomic weight missing')
            catalogue.append(dict(Species=name, Formula=formula, Phase='molecular', Ready=ready,
                                  Reason=reason, NASA=mapped))
            if ready:
                entry = dict(formula=formula, nasa=mapped, elements=elements)
                self.molecular[name] = entry
                self.organic[name] = entry

        # Every NASA record stays usable as a dry-feed species; solid records are pure solids.
        covered = {entry.get('nasa') for entry in self.molecular.values()}
        for name, data in nasa_db.items():
            elements = data.get('elements')
            if not elements or not all(e in AW for e in elements):
                continue
            key = name if name not in original_molecules else name.casefold()
            entry = dict(formula=name, nasa=name, elements=elements)
            if name not in covered:
                self.molecular.setdefault(key, entry)
                catalogue.append(dict(Species=key, Formula=name, Phase='molecular', Ready=False,
                                      Reason='NASA available; no resolved SAFT mapping. '
                                             'Supply species_overrides or EOS data.', NASA=name))
            if data.get('phase') == 'S':
                self.solid.setdefault(key, entry)
                catalogue.append(dict(Species=key, Formula=name, Phase='solid', Ready=True,
                                      Reason='NASA solid', NASA=name))
        for name, entry in self.solid.items():
            catalogue.append(dict(Species=name, Formula=entry['formula'], Phase='solid',
                                  Ready=entry['nasa'] in nasa_db, Reason='NASA solid',
                                  NASA=entry['nasa']))

        # Aqueous ions: charged ePC-SAFT rows joined to SUPCRT records by charge and identity.
        xml = list(databases.supcrt.iter('Species'))
        ready_epc, ready_hkf = [], []
        charged = epcsaft[epcsaft['record_type'].astype(str).eq('species_parameter')
                          & epcsaft['Charge'].fillna(0).ne(0)]
        for name, rows in charged.groupby('Species', sort=False):
            name = str(name).strip()
            charge = int(rows.Charge.dropna().iloc[0])
            formula = str(rows.Formula.dropna().iloc[0])
            self.ion_charge[name] = charge
            self.ion_saft_formula[name] = formula
            matches = [node for node in xml if node.findtext('Type') == 'Aqueous'
                       and float(node.findtext('Charge') or 0) == charge
                       and any(_ion_identity(candidate, charge) == _ion_identity(name, charge)
                               for candidate in (node.findtext('Name'), node.findtext('Formula')))]
            self.supcrt_ion_aliases[name] = [node.findtext('Name') for node in matches]
            epc_data = first_complete_ion_epcsaft_set(epcsaft, name, formula, charge)
            hkf_data = (first_complete_hkf_set(databases.supcrt, name, self.supcrt_ion_aliases[name])
                        if matches else None)
            node = next((node for node in matches if read_hkf_candidate(node)), None)
            elements = ({e: float(v) for e, v in re.findall(ELEMENT_PATTERN,
                                                             node.findtext('Elements') or '')}
                        if node is not None else None)
            ready = bool(epc_data and hkf_data and elements and all(e in AW for e in elements))
            catalogue.append(dict(Species=name, Formula=formula, Phase='aqueous ion', Ready=ready,
                                  Reason='ready' if ready else
                                  'Complete matching ionic EOS/HKF/element data missing', NASA=None))
            if ready:
                self.ion[name] = dict(formula=formula, charge=charge, supcrt=hkf_data['SUPCRT name'],
                                      elements=elements)
                ready_epc.append(epc_data)
                ready_hkf.append(hkf_data)
        self.ion_epcsaft_table = pd.DataFrame(ready_epc)
        self.ion_hkf_table = pd.DataFrame(ready_hkf)
        # SUPCRT-only records are listed so users can see what exists.
        ion_names = {row['Species'] for row in catalogue if row['Phase'] == 'aqueous ion'}
        for node in xml:
            name, kind = node.findtext('Name'), node.findtext('Type')
            if name in ion_names:
                continue
            catalogue.append(dict(Species=name, Formula=node.findtext('Formula'),
                                  Phase=f'SUPCRT {kind}', Ready=False,
                                  Reason='SUPCRT reference present; matching EOS/reference '
                                         'adapter required', NASA=None))
        self.catalogue = pd.DataFrame(catalogue).drop_duplicates(['Species', 'Phase'])
        self._register_supcrt_minerals(databases)

    def _register_supcrt_minerals(self, databases):
        """Single-range Maier-Kelley minerals from supcrt92.xml become pure solids.

        Polymorphic records (phase transitions) are listed as unavailable rather
        than extrapolated. Mineral <Thermo/HKF> records are not aqueous HKF records.
        """
        rows = []
        for node in databases.supcrt.iter('Species'):
            if node.findtext('Type') != 'Mineral':
                continue
            name = node.findtext('Name')
            hkf = node.find('./Thermo/HKF')
            reason, params = '', {}
            try:
                if hkf is None:
                    raise ValueError('Mineral HKF data missing')
                if int(hkf.findtext('NumPhaseTrans') or 0):
                    raise ValueError('Polymorphic mineral: transition adapter not implemented')
                for tag in ('Gf', 'Hf', 'Sr', 'Vr', 'a', 'b', 'c', 'Tmax'):
                    params[tag] = float(hkf.findtext(tag))
                if (not np.all(np.isfinite(list(params.values()))) or params['Vr'] <= 0
                        or params['Tmax'] < 298.15):
                    raise ValueError('Invalid/missing mineral parameters')
                for tag in ('Gf', 'Hf', 'Sr', 'a', 'b', 'c'):
                    if 'cal' not in hkf.find(tag).get('units', ''):
                        raise ValueError(f'Unexpected units for mineral {tag}')
                if hkf.find('Vr').get('units') != 'cm3/mol':
                    raise ValueError('Expected Vr in cm3/mol (do not use mislabeled MolarVolume)')
                elements = {e: float(v) for e, v in re.findall(ELEMENT_PATTERN,
                                                               node.findtext('Elements') or '')}
                if not elements or any(e not in ATOMIC_WEIGHTS for e in elements):
                    raise ValueError('Element composition or atomic weight missing')
                formula = node.findtext('Formula') or ''
                self.solid[name] = dict(
                    formula=node.findtext('Formula'), elements=elements, adapter='supcrt_mineral',
                    nasa=None, supcrt_mineral=params, source='supcrt92.xml',
                    references=hkf.findtext('References'), date=hkf.findtext('Date'),
                    hydrate=('*' in formula and 'H2O' in formula))
            except (TypeError, ValueError) as exc:
                reason = str(exc)
            rows.append(dict(Species=name, Formula=node.findtext('Formula'), Phase='solid',
                             Ready=not bool(reason),
                             Reason=reason or 'SUPCRT mineral: single-range Maier-Kelley', NASA=None))
        names = {r['Species'] for r in rows}
        keep = ~(self.catalogue.Species.isin(names)
                 & self.catalogue.Phase.isin(['SUPCRT Mineral', 'solid']))
        self.catalogue = pd.concat([self.catalogue[keep], pd.DataFrame(rows)], ignore_index=True)

    def copy(self):
        return copy.deepcopy(self)
