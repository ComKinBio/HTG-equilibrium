"""Write test/data/python_results/: a run saved by the Python app, for the Julia tests.

    uv run python julia/HTGEquilibrium/test/dump_results.py

The run is made of hand-made Results (no solving), so the files are small and
cover what the tables can hold: ions with molality, an empty phase, a failed
point whose message has a comma, quotes and a line break, a point with the
aqueous phase removed and a dropped carbon bound, missing descent_RT, and
points with two and three results (alternatives.csv). Written:

    spec.json          the Results (Julia builds the same ones and must write the same files)
    run/               SweepResult.save + write_choices (summary.csv, species.csv, phases.csv,
                       alternatives.csv, case.toml, choices.toml)
    python_read.json   what Python reads back: read_choices, and htg.gui.runs.load_saved's
                       tables (choices applied, and as saved) with their dtypes
"""
import json
import math
import shutil
from pathlib import Path

import numpy as np

import htg
from htg.gui.runs import load_saved
from htg.results import Alternative, Result, SweepResult, read_choices, write_choices

HERE = Path(__file__).resolve().parent
OUT = HERE / 'data' / 'python_results'
TEMPLATE = Path(htg.__file__).parent / 'templates' / 'case_template.toml'

AW = {'C': 12.011, 'H': 1.008, 'O': 15.999, 'N': 14.007, 'Na': 22.99}
SPECIES = {  # (phase, species): (formula, elements, charge)
    ('gas', 'water'): ('H2O', {'H': 2, 'O': 1}, 0),
    ('gas', 'carbon dioxide'): ('CO2', {'C': 1, 'O': 2}, 0),
    ('gas', 'methane'): ('CH4', {'C': 1, 'H': 4}, 0),
    ('gas', 'hydrogen'): ('H2', {'H': 2}, 0),
    ('aqueous', 'water'): ('H2O', {'H': 2, 'O': 1}, 0),
    ('aqueous', 'carbon dioxide'): ('CO2', {'C': 1, 'O': 2}, 0),
    ('aqueous', 'Na+'): ('Na+', {'Na': 1}, 1),
    ('aqueous', 'HCO3-'): ('HCO3-', {'H': 1, 'C': 1, 'O': 3}, -1),
    ('organic_liquid', 'benzene'): ('C6H6', {'C': 6, 'H': 6}, 0),
    ('solid', 'carbon'): ('C', {'C': 1}, 0),
}
ELEMENTS = ['C', 'H', 'N', 'Na', 'O']


def _amounts(scale, aqueous=True):
    a = {('gas', 'water'): 0.4 * scale, ('gas', 'carbon dioxide'): 0.5, ('gas', 'methane'): 1 / 3,
         ('gas', 'hydrogen'): 1e-20, ('aqueous', 'water'): 12.5 if aqueous else 0.0,
         ('aqueous', 'carbon dioxide'): 0.05 * scale if aqueous else 0.0,
         ('aqueous', 'Na+'): 0.01 if aqueous else 0.0, ('aqueous', 'HCO3-'): 0.01 if aqueous else 0.0,
         ('organic_liquid', 'benzene'): 0.0, ('solid', 'carbon'): 0.1 / scale}
    return a


def _phase_gibbs(amounts, offset):
    out = {}
    for (phase, _), v in amounts.items():
        out[phase] = out.get(phase, 0.0) - 1000.0 * v - offset
    return out


def make_result(T, aqueous=True, alternatives=(), descent=True, dropped=()):
    amounts = _amounts(T / 300.0, aqueous)
    gibbs = -5e6 - 1234.5678 * T
    feed_C = 1.0
    gas_C = sum(v * SPECIES[p][1].get('C', 0) for p, v in amounts.items() if p[0] == 'gas')
    optimizer = {'method': 'ipopt', 'nfev': 40}
    if descent:
        optimizer['descent_RT'] = 1.5e-05 * T / 300
    alts = [Alternative(label, gibbs + dG, _amounts(scale, aqueous), n, starts, _phase_gibbs(amounts, dG))
            for label, dG, scale, n, starts in alternatives]
    return Result(
        temperature_C=float(T), pressure_bar=240.0, amounts=amounts, initial_amounts=dict(amounts),
        gibbs_J_per_mol_C=gibbs, phase_gibbs_J_per_mol_C=_phase_gibbs(amounts, 0.0),
        balance_residual=3.552713678800501e-15, charge_residual=0.0, iterations=23 + int(T) // 50,
        eos_trial_errors=[], feed_C_mol_per_g_dry=0.0368, elements=ELEMENTS,
        feed_atoms=np.array([feed_C, 2.0, 0.0, 0.01, 1.5]), product_atoms=np.array([feed_C, 2.0, 0.0, 0.01, 1.5]),
        gas_carbon_fraction=gas_C / feed_C, carbon_bound_violation=1.1102230246251565e-16,
        phase_rules={'aqueous_disabled': not aqueous, 'dropped_constraints': list(dropped)},
        optimizer=optimizer, atomic_weights=AW,
        species_data={p: {'formula': f, 'elements': e, 'charge': c} for p, (f, e, c) in SPECIES.items()},
        alternatives=alts)


ERROR = 'optimization / final validation: RuntimeError: Iteration limit, "x" reached\nsecond line'


def spec_of(T, P, r):
    if isinstance(r, str):
        return {'T_C': T, 'P_bar': P, 'error': r}
    pairs = lambda d: [[p, s, v] for (p, s), v in d.items()]  # noqa: E731
    return {'T_C': T, 'P_bar': P, 'gibbs_J_per_mol_C': r.gibbs_J_per_mol_C,
            'amounts': pairs(r.amounts), 'phase_gibbs_J_per_mol_C': r.phase_gibbs_J_per_mol_C,
            'balance_residual': r.balance_residual, 'charge_residual': r.charge_residual,
            'iterations': r.iterations, 'feed_C_mol_per_g_dry': r.feed_C_mol_per_g_dry, 'elements': r.elements,
            'feed_atoms': list(r.feed_atoms), 'gas_carbon_fraction': r.gas_carbon_fraction,
            'carbon_bound_violation': r.carbon_bound_violation, 'phase_rules': r.phase_rules,
            'optimizer': r.optimizer, 'atomic_weights': r.atomic_weights,
            'species_data': [[p, s, d] for (p, s), d in r.species_data.items()],
            'alternatives': [{'label': a.label, 'gibbs': a.gibbs, 'amounts': pairs(a.amounts), 'n_starts': a.n_starts,
                              'start_labels': a.start_labels, 'phase_gibbs': a.phase_gibbs}
                             for a in r.alternatives]}


def _plain(value):
    if isinstance(value, (float, np.floating)):
        return None if math.isnan(value) else float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def table_json(frame):
    return {'columns': list(frame.columns), 'dtypes': [str(t) for t in frame.dtypes],
            'rows': [[_plain(v) for v in row] for row in frame.itertuples(index=False)]}


def main():
    case = htg.load_case(TEMPLATE).replace(temperature_C=[300.0, 350.0, 400.0, 450.0], pressure_bar=240.0,
                                           aqueous_cutoff_C=None)
    points = [
        (300.0, 240.0, make_result(300)),
        (350.0, 240.0, ERROR),
        (400.0, 240.0, make_result(400, aqueous=False, descent=False, dropped=['aqueous carbon'],
                                   alternatives=[('A', 0.0, 400 / 300, 3, ['lp', 'random 1', 'random 2']),
                                                 ('B', 120.0, 1.5, 2, ['phase gas', 'random 3']),
                                                 ('C', 900.25, 0.5, 1, ['random 4'])])),
        (450.0, 240.0, make_result(450, alternatives=[('A', 0.0, 1.5, 1, ['lp']),
                                                      ('B', 30.0, 2.0, 1, ['random 1'])])),
    ]
    if OUT.exists():
        shutil.rmtree(OUT)
    run = OUT / 'run'
    SweepResult(case, points).save(run)
    write_choices(run, {(400.0, 240.0): 'B', (450, 240): 'A', (300.0, 240.0): 'Z'})
    saved = load_saved(run, OUT)
    read = {'choices': [[T, P, label] for (T, P), label in read_choices(run).items()],
            'label': saved.label, 'flagged_points': saved.flagged_points,
            'applied': {name: table_json(getattr(saved, name)) for name in ('summary', 'species', 'phases')},
            'base': {name: table_json(t) for name, t in zip(('summary', 'species', 'phases'), saved.base)},
            'alternatives': table_json(saved.alternatives)}
    (OUT / 'python_read.json').write_text(json.dumps(read, separators=(',', ':')))
    spec = {'case_temperature_C': case.temperature_C, 'points': [spec_of(*p) for p in points]}
    (OUT / 'spec.json').write_text(json.dumps(spec, separators=(',', ':')))
    size = sum(f.stat().st_size for f in OUT.rglob('*') if f.is_file())
    print(f'{OUT}: {size / 1000:.0f} kB')


if __name__ == '__main__':
    main()
