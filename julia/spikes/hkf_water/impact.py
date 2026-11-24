"""What option 2 does to the app's results: rerun the electrolyte reference
checks with the option-2 water properties and ion G patched in.

    uv run python julia/spikes/hkf_water/impact.py

`standard_states.supcrt_water` and `standard_states.hkf_ion_g0` are replaced by
lookups into julia_results.json (IAPWS-95 density + Julia JN91/Born/HKF), at the
(T, P) points of the comparison grid. Everything else (ePC-SAFT, optimizer) is
the current Python app. Compares against tests/reference/data with the
tolerances of tests/test_reference_thermo.py and tests/test_reference_solves.py.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'tests' / 'reference'))

import htg  # noqa: E402
from htg import standard_states as ss  # noqa: E402
from htg.constants import ATOMIC_WEIGHTS  # noqa: E402
from htg.phases import ElectrolytePhase  # noqa: E402
from convert import case_from_legacy, load, pairs  # noqa: E402

JL = json.loads((HERE / 'julia_results.json').read_text())
TABLE = {(round(p['T_C'], 3), round(p['P_bar'], 3)): p for p in JL['points']}
ORIGINAL = {'water': ss.supcrt_water, 'hkf': ss.hkf_ion_g0}


def key(T, P):
    return round(float(T) - 273.15, 3), round(float(P) / 1e5, 3)


def patched_water(T, P):
    return dict(TABLE[key(T, P)]['water_iapws'])


def patched_hkf(hkf_table, ions, T, P):
    mu = TABLE[key(T, P)]['mu0_ion']
    return np.array([mu[i]['iapws'] for i in ions])


def patch(on):
    ss.supcrt_water = patched_water if on else ORIGINAL['water']
    ss.hkf_ion_g0 = patched_hkf if on else ORIGINAL['hkf']


def phase_level(model):
    print('## Electrolyte phase at the thermo.json states (option 2 vs reference)\n')
    print('| state | selection | max rel d mu0_ion | max abs d lnphi_inf | max abs d ln phi (points) | max rel d G (points) |')
    print('|---|---|---|---|---|---|')
    reg = model.registry
    for state in load('thermo')['states']:
        T, P = state['T_C'] + 273.15, state['P_bar'] * 1e5
        for label, ref in state.get('electrolyte', {}).items():
            sel = ref['selection']
            neutrals = list(dict.fromkeys(reg.aqueous_reference_neutrals + [s for s in sel if s not in reg.ion]))
            resolved = reg.resolve_parameters(model.databases, T, neutrals, ATOMIC_WEIGHTS)
            patch(True)
            try:
                ph = ElectrolytePhase(model, sel, T, P, resolved, ATOMIC_WEIGHTS)
            finally:
                patch(False)
            dmu = np.max(np.abs(np.array(ph.mu0_ion) / np.array(ref['mu0_ion']) - 1))
            dlinf = np.max(np.abs(np.array(ph.lnphi_infinity) - np.array(ref['lnphi_infinity'])))
            dG, dlnphi = 0.0, 0.0
            for point in ref['points']:
                if 'error' in point:
                    continue
                full = np.array([point['amounts'].get(s, 0.0) for s in ph.all_species])
                dG = max(dG, abs(ph.gibbs_full(full) / point['G'] - 1))
                dlnphi = max(dlnphi, np.max(np.abs(np.log(ph.phi(full)) - np.array(point['ln_phi']))))
            print(f"| {state['T_C']:g} C, {state['P_bar']:g} bar | {label} | {dmu:.1e} | {dlinf:.1e} | "
                  f"{dlnphi:.1e} | {dG:.1e} |")


CASES = ['fig78_300C_ipopt', 'full_text_250C', 'full_text_gas_min_300C', 'configurable_300C_ipopt',
         'fig78_380C']


def solve_level(models):
    print('\n## Reference solves with an aqueous electrolyte phase\n')
    print('| case | T C | variant | rel d G (test: 1e-9) | max abs d n (mol/mol C) | '
          'max rel d n, n > 1e-4 | d gas C fraction (test: 1e-7) |')
    print('|---|---|---|---|---|---|---|')
    ref_cases = {c['name']: c for c in load('solves')['cases']}
    notes = []
    for name in CASES:
        entry = ref_cases[name]
        ref = entry['result']
        case = case_from_legacy(entry)
        if name == 'fig78_380C':   # as tests/test_reference_solves.py (FD_STOPS_EARLY)
            case = case.replace(solver='ipopt')
        expected = pairs(ref['amounts'])
        for variant in ('pychnosz', 'option 2'):
            patch(variant == 'option 2')
            try:
                r = models[entry['model']].solve(case)
            finally:
                patch(False)
            n = np.array([r.amounts.get(k, 0.0) for k in expected])
            e = np.array(list(expected.values()))
            big = np.abs(e) > 1e-4
            print(f"| {name} | {entry['settings']['temperature_C']:g} | {variant} | "
                  f"{r.gibbs_J_per_mol_C / ref['gibbs_J_per_mol_feed_C'] - 1:.1e} | "
                  f"{np.max(np.abs(n - e)):.1e} | {np.max(np.abs(n[big] / e[big] - 1)):.1e} | "
                  f"{r.gas_carbon_fraction - ref['gas_carbon_fraction_actual']:.1e} |")
            if variant == 'option 2':
                rel = np.where(big, np.abs(n / np.where(big, e, 1) - 1), 0)
                keys = list(expected)
                notes.append(f'{name}: ' + ', '.join(
                    f'{keys[i][1]} ({keys[i][0]}) {e[i]:.4g} -> {n[i]:.4g}'
                    for i in np.argsort(-rel)[:4] if rel[i] > 1e-5))
    print('\nLargest relative amount changes with option 2 (n > 1e-4 mol/mol C):\n')
    for line in notes:
        print('- ' + line)


def main():
    frozen = ROOT / 'tests' / 'reference' / 'databases'
    models = {'base': htg.Model.load(frozen, nasa_names={}), 'acetone': htg.Model.load(frozen)}
    phase_level(models['base'])
    solve_level(models)


if __name__ == '__main__':
    main()
