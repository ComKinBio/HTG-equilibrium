"""Generate the reference results from the FROZEN legacy model.

    uv run python tests/reference/generate.py [--data DIR] [--only thermo|catalogue|solves]

Writes tests/reference/data/*.json. Run this only to extend the reference set
with new cases; the frozen model in legacy/ must never be edited, so the stored
results always describe the original notebook implementation.
"""
import argparse
import importlib.util
import json
import sys
import time
import traceback
import types
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
# The legacy notebooks ran from the repository root, where juliapkg.json (Ipopt) lives.
sys.path.insert(0, str(HERE.parents[1]))
import legacy_cases  # noqa: E402

OUT = HERE / 'data'


def load_legacy(data_dir, species_overrides=None):
    spec = importlib.util.spec_from_file_location('htg_optimizer', HERE / 'legacy' / 'htg_optimizer.py')
    optimizer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(optimizer)
    sys.modules['htg_optimizer'] = optimizer  # the legacy model does `import htg_optimizer`
    path = HERE / 'legacy' / 'htg_model_legacy.py'
    module = types.ModuleType('htg_model_legacy')
    exec(compile(path.read_text(encoding='utf-8'), str(path), 'exec'), module.__dict__)
    module.initialize_model(database_paths(data_dir), species_overrides=species_overrides)
    return module


def database_paths(data_dir):
    return {'nasa': data_dir / 'Thermdat_WasteMech.txt', 'pcsaft': data_dir / 'PROJECT_PCSAFT.csv',
            'epcsaft': data_dir / 'PROJECT_EPCSAFT.csv', 'supcrt': data_dir / 'supcrt92.xml'}


def jsonable(value):
    if isinstance(value, dict):
        if value and all(isinstance(k, tuple) for k in value):
            return [[*k, jsonable(v)] for k, v in value.items()]
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return [jsonable(v) for v in value.tolist()]
    if isinstance(value, (np.floating, float)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def error_text(exc):
    return f'{type(exc).__name__}: {exc}'


def write(name, payload):
    OUT.mkdir(exist_ok=True)
    (OUT / f'{name}.json').write_text(json.dumps(jsonable(payload), indent=1, sort_keys=False) + '\n')
    print(f'wrote {name}.json', flush=True)


# ---------------------------------------------------------------- catalogue
def catalogue(m):
    cat = m.available_species(include_unavailable=True)
    registries = {}
    for label, reg in [('molecular', m.MOLECULAR_SPECIES), ('organic', m.ORGANIC_SPECIES),
                       ('solid', m.SOLID_DATA), ('ion', m.ION_DATA)]:
        registries[label] = {name: {k: d.get(k) for k in ('formula', 'nasa', 'elements', 'charge',
                                                          'supcrt', 'adapter', 'supcrt_mineral',
                                                          'hydrate')
                                    if k in d}
                             for name, d in reg.items()}
    return {
        'catalogue': cat.astype(object).where(cat.notna(), None).to_dict('records'),
        'registries': registries,
        'ion_epcsaft_table': m.ion_epcsaft_table.astype(object).where(m.ion_epcsaft_table.notna(), None).to_dict('records'),
        'ion_hkf_table': m.ion_hkf_table.astype(object).where(m.ion_hkf_table.notna(), None).to_dict('records'),
        'mineral_inventory': m.mineral_inventory().to_dict('records'),
        'ion_product_candidates_400C': m.ion_product_candidates(temperature_C=400.0).to_dict('records'),
    }


# ---------------------------------------------------------------- thermodynamics
PARAM_KEYS = ('Formula', 'Source', 'Database row', 'Model', 'RefID', 'm', 'sigma', 'epsilon',
              'epsilonAB', 'kappaAB', 'association_scheme', 'Mw')


def thermo(m):
    rng = np.random.default_rng(20261110)
    ready = m.available_species()
    ready_molecular = ready.loc[ready.Phase.eq('molecular'), 'Species'].tolist()
    states = []
    for T_C, P_bar, electrolyte in legacy_cases.STATES:
        t0 = time.perf_counter()
        T, P = T_C + 273.15, P_bar * 1e5
        m.T, m.P = T, P
        state = {'T_C': T_C, 'P_bar': P_bar}
        g = {}
        for name in sorted(m.ALL_MOLECULAR_DATA):
            try:
                g[name] = m.nasa7_properties(m.nasa_db, m.get_nasa_name(name), T)[3]
            except Exception as exc:
                g[name] = error_text(exc)
        state['molecular_g0'] = g
        solids = {}
        for name in sorted(m.SOLID_DATA):
            try:
                solids[name] = float(m.solid_standard_chemical_potentials([name])[0][0])
            except Exception as exc:
                solids[name] = error_text(exc)
        state['solid_g0'] = solids
        params = {}
        for name in ready_molecular:
            try:
                d = m.resolve_molecular_parameters(T, [name])[name]
                params[name] = {k: d.get(k) for k in PARAM_KEYS}
            except Exception as exc:
                params[name] = error_text(exc)
        state['saft_parameters'] = params

        fluids = {}
        for label, (names, root) in legacy_cases.FLUID_PHASES.items():
            m._refresh_state(T_C, P_bar, electrolyte=False, phase_species={'gas': names})
            mu0, _ = m.nasa_standard_chemical_potentials(names)
            model, _ = m.build_pcsaft_model(names)
            points = []
            water = (1e-3, 0.5) if root == 'vapour' else (2.0, 20.0)
            for n in legacy_cases.random_compositions(names, rng, water=water):
                vector = np.array([n[s] for s in names])
                point = {'amounts': n}
                try:
                    point['G'] = m.molecular_gibbs(vector, names, mu0, root)
                    phis = m.cl.fugacity_coefficient(model, P, T, np.maximum(vector, 1e-14), phase=root)
                    point['ln_phi'] = np.log(np.asarray(phis, dtype=float)).tolist()
                except Exception as exc:
                    point['error'] = error_text(exc)
                points.append(point)
            fluids[label] = {'species': names, 'root': root, 'mu0': mu0, 'points': points}
        state['fluid_phases'] = fluids

        if electrolyte:
            electro = {}
            for label, selection in legacy_cases.ELECTROLYTE_SELECTIONS.items():
                m.configure_aqueous_species({'aqueous': selection})
                m._refresh_state(T_C, P_bar, electrolyte=True, phase_species={'aqueous': selection})
                charges = {s: m.ION_DATA[s]['charge'] for s in selection if s in m.ION_DATA}
                points = []
                for n in legacy_cases.random_compositions(selection, rng, ions=charges):
                    full = np.array([n.get(s, 0.0) for s in m.AQUEOUS_SPECIES])
                    point = {'amounts': n}
                    try:
                        point['G'] = m.aqueous_gibbs(full)
                        point['ln_phi'] = np.log(m.phi_electrolyte(full)).tolist()
                    except Exception as exc:
                        point['error'] = error_text(exc)
                    points.append(point)
                electro[label] = {
                    'selection': selection, 'aqueous_species': list(m.AQUEOUS_SPECIES),
                    'water': m.water_operating, 'mu0_ion': m.mu0_ion,
                    'mu0_neutral': m.mu0_aq_neutral, 'lnphi_infinity': m.lnphi_infinity,
                    'points': points,
                }
            m.configure_aqueous_species({'aqueous': m._DEFAULT_AQ_NEUTRALS + m._DEFAULT_AQ_IONS})
            state['electrolyte'] = electro
        states.append(state)
        print(f'state {T_C:g} C / {P_bar:g} bar done in {time.perf_counter() - t0:.0f} s', flush=True)
    return {'states': states}


# ---------------------------------------------------------------- solves
RESULT_KEYS = ('temperature_C', 'pressure_bar', 'amounts', 'initial_amounts',
               'gibbs_J_per_mol_feed_C', 'phase_gibbs_J_per_mol_feed_C', 'balance_residual',
               'charge_residual', 'iterations', 'feed_C_mol_per_g_dry', 'elements', 'feed_atoms',
               'product_atoms', 'gas_carbon_fraction_actual', 'gas_carbon_min_violation',
               'phase_allocation', 'eos_trial_errors')


def solves(models, data_dir):
    cases = legacy_cases.build_cases(models['base'], data_dir)
    results, finished = [], {}
    for case in cases:
        m = models[case['model']]
        settings = dict(case['settings'])
        if case['continue_from']:
            settings['initial_phase_amounts'] = dict(finished[case['continue_from']]['amounts'])
        saved_AW = m.AW
        t0 = time.perf_counter()
        entry = {'name': case['name'], 'model': case['model'],
                 'atomic_weights': case['atomic_weights'], 'settings': settings,
                 'continue_from': case['continue_from']}
        try:
            m.AW = {**m.AW, **case['atomic_weights']}
            r = m.solve_equilibrium(settings)
            finished[case['name']] = r
            entry['result'] = {k: r[k] for k in RESULT_KEYS}
            if 'optimizer' in r:
                o = r['optimizer']
                entry['result']['optimizer'] = {k: o[k] for k in ('method', 'starts', 'descent_RT',
                                                                  'complementarity_J',
                                                                  'stationarity_major_RT',
                                                                  'tangent_plane_RT') if k in o}
        except Exception as exc:
            entry['error'] = error_text(exc)
            traceback.print_exc()
        finally:
            m.AW = saved_AW
        entry['seconds'] = time.perf_counter() - t0
        status = 'error' if 'error' in entry else f"G={entry['result']['gibbs_J_per_mol_feed_C']:.6f}"
        print(f"{case['name']}: {status} ({entry['seconds']:.0f} s)", flush=True)
        results.append(entry)
    return {'cases': results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=HERE.parents[1] / 'data')
    parser.add_argument('--only', choices=['catalogue', 'thermo', 'solves'])
    args = parser.parse_args()
    data_dir = args.data.resolve()
    t0 = time.perf_counter()
    base = load_legacy(data_dir)
    print(f'legacy model initialized in {time.perf_counter() - t0:.0f} s', flush=True)
    acetone = load_legacy(data_dir, species_overrides={'acetone': {'nasa': 'CH3COCH3'}})
    if args.only in (None, 'catalogue'):
        write('catalogue', {'base': catalogue(base), 'acetone': catalogue(acetone)})
    # Solves first, so they start from a freshly initialized model as in the notebooks.
    if args.only in (None, 'solves'):
        write('solves', solves({'base': base, 'acetone': acetone}, data_dir))
    if args.only in (None, 'thermo'):
        write('thermo', thermo(base))


if __name__ == '__main__':
    main()
