"""Write the user species fixture for test/test_registry.jl.

    uv run python julia/HTGEquilibrium/test/dump_user_species.py

1. Builds a user species folder, test/data/user_species/, with the Python app's
   UserDatabase (the specs of tests/test_user_species.py: a PC-SAFT molecule,
   an ePC-SAFT molecule, an override of methane, an ion copied from SUPCRT92,
   a Maier-Kelley mineral, a NASA solid, an incomplete species, a hidden
   shipped species).
2. Loads the model from tests/reference/databases with that folder and writes
   what the Julia port must reproduce to test/data/user_species_python.json.
"""
import json
import math
import shutil
from pathlib import Path

import numpy as np

import htg
from htg import standard_states as ss
from htg.user_species import SpeciesSpec, UserDatabase, spec_from_model

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
FROZEN = REPO / 'tests' / 'reference' / 'databases'
FOLDER = HERE / 'data' / 'user_species'
OUT = HERE / 'data' / 'user_species_python.json'
T, P = 623.15, 240e5
PARAM_KEYS = ('Species', 'Formula', 'Source', 'Database row', 'Model', 'RefID', 'm', 'sigma', 'epsilon',
              'epsilonAB', 'kappaAB', 'association_scheme', 'Mw')


def plain(value):
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [plain(v) for v in value]
    if isinstance(value, (np.floating, float)):
        return None if math.isnan(value) else float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, set):
        return sorted(value)
    return value


def records(frame):
    return plain(frame.astype(object).where(frame.notna(), None).to_dict('records'))


def nasa_coefficients(model, record):
    r = model.databases.nasa[record]
    return {'T_low': r['T_low'], 'T_mid': r['T_mid'], 'T_high': r['T_high'],
            'high': list(r['high']), 'low': list(r['low'])}


def build(shipped):
    if FOLDER.exists():
        shutil.rmtree(FOLDER)
    db = UserDatabase.load(FOLDER)
    db.save_species(SpeciesSpec(
        kind='molecule', name='Testanol', formula='C9H12O3', thermo='nasa',
        nasa=nasa_coefficients(shipped, 'CH4'), model='PC-SAFT', m=3.1, sigma=3.6, epsilon=250.0,
        association_scheme='2B', epsilonAB=2000.0, kappaAB=0.02, kij=[('methane', 0.01)],
        authors='A. Author', year='2024', doi='10.1/x'), shipped)
    db.save_species(SpeciesSpec(
        kind='molecule', name='Testol', formula='C5H12O2', thermo='nasa',
        nasa=nasa_coefficients(shipped, 'CH3OH'), model='ePC-SAFT', m=2.5, sigma=3.4, epsilon=260.0,
        association_scheme='2B', epsilonAB=2500.0, kappaAB=0.03), shipped)
    methane = spec_from_model(shipped, 'methane')
    methane.m = methane.m * 1.5
    db.save_species(methane, shipped)
    db.save_species(SpeciesSpec(kind='ion', name='Ag+', formula='Ag+', charge=1, thermo='hkf_record',
                                supcrt_record='Ag+', m=1.0, sigma=2.9, epsilon=200.0,
                                dielectric_constant=8.0), shipped)
    db.save_species(SpeciesSpec(kind='solid', name='Testite', formula='CaCO3', thermo='maier_kelley',
                                maier_kelley={'Gf': -269880.0, 'Hf': -288552.0, 'Sr': 22.15, 'Vr': 36.934,
                                              'a': 24.98, 'b': 0.00524, 'c': -620000.0, 'Tmax': 1200.0}),
                    shipped)
    db.save_species(SpeciesSpec(kind='solid', name='mychar', formula='C', thermo='nasa',
                                nasa=nasa_coefficients(shipped, 'CHAR')), shipped)
    db.save_species(SpeciesSpec(kind='molecule', name='halfdone', formula='C2H6O', thermo='nasa',
                                nasa={'T_low': 300, 'T_mid': 1000}, m=None, sigma=3.0, epsilon=None),
                    shipped)
    UserDatabase.load(FOLDER).hide_species('ethane')


def main():
    shipped = htg.Model.load(FROZEN)
    build(shipped)
    model = htg.Model.load(FROZEN, user_data=FOLDER)
    reg = model.registry
    touched = set(model.user_info) | {'ethane', 'testanol', 'testol'}
    catalogue = records(model.species(include_unavailable=True))
    molecules = ['testanol', 'testol', 'methane']
    out = {
        'T': T, 'P': P,
        'user_info': plain(model.user_info),
        'catalogue_length': len(catalogue),
        'catalogue_rows': [r for r in catalogue if r['Species'] in touched],
        'registry_names': {k: list(getattr(reg, k)) for k in ('molecular', 'organic', 'solid', 'ion')},
        'entries': {name: plain(reg.any_entry(name)) for name in
                    ['testanol', 'testol', 'methane', 'Ag+', 'Testite', 'mychar']},
        'saft_parameters': {name: plain({k: d.get(k) for k in PARAM_KEYS}) for name, d in
                            reg.resolve_parameters(model.databases, T, molecules).items()},
        'binary_interactions': plain(reg.binary_interactions),
        'ion_epcsaft_table': records(reg.ion_epcsaft_table),
        'ion_hkf_table': records(reg.ion_hkf_table),
        'nasa_g0': dict(zip(molecules, plain(ss.nasa_g0(reg, model.databases.nasa, molecules, T)))),
        'solid_g0': dict(zip(['Testite', 'mychar'],
                             plain(ss.solid_g0(reg, model.databases.nasa, ['Testite', 'mychar'], T, P)))),
        'hkf_ion_g0': dict(zip(['Ag+', 'Na+'],
                               plain(ss.hkf_ion_g0(reg.ion_hkf_table, ['Ag+', 'Na+'], 573.15, 240e5)))),
    }
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False) + '\n', encoding='utf-8')
    print(f'wrote {FOLDER.relative_to(REPO)} and {OUT.relative_to(REPO)}')


if __name__ == '__main__':
    main()
