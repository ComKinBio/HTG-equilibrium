"""Write test/data/solves_python_ipopt.json: the Python app's IPOPT solution of every
solved point of tests/reference/data/solves.json.

Most reference points were solved with the original finite-difference SLSQP
("slsqp-fd"), which the Julia version does not have; the Julia tests solve them
with IPOPT and compare with this file as well (Python's IPOPT used a
finite-difference Hessian, Julia's is exact; the converged points must agree).

    uv run python julia/HTGEquilibrium/test/dump_solves.py
"""
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'tests' / 'reference'))

import htg  # noqa: E402
from convert import case_from_legacy, load  # noqa: E402

FROZEN = ROOT / 'tests' / 'reference' / 'databases'
OUT = Path(__file__).resolve().parent / 'data' / 'solves_python_ipopt.json'


def main():
    models = {'base': htg.Model.load(FROZEN, nasa_names={}), 'acetone': htg.Model.load(FROZEN)}
    out = {}
    for entry in load('solves')['cases']:
        if 'result' not in entry:
            continue
        case = case_from_legacy(entry)
        if case.solver == 'slsqp-fd':
            case = case.replace(solver='ipopt')
        t0 = time.perf_counter()
        r = models[entry['model']].solve(case)
        o = r.optimizer
        out[entry['name']] = {
            'solver': case.solver,
            'gibbs_J_per_mol_C': r.gibbs_J_per_mol_C,
            'phase_gibbs_J_per_mol_C': r.phase_gibbs_J_per_mol_C,
            'amounts': [[p, s, v] for (p, s), v in r.amounts.items()],
            'iterations': r.iterations,
            'starts': [{k: run.get(k) for k in ('start', 'G', 'message', 'nit')} for run in o['starts']],
            'stationarity_major_RT': o.get('stationarity_major_RT'),
            'polish_used': o.get('polish', {}).get('used'),
            'seconds': time.perf_counter() - t0,
        }
        print(f"{entry['name']}: G = {r.gibbs_J_per_mol_C!r}", flush=True)
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(out, indent=0) + '\n')
    print(f'{OUT}: {len(out)} points, {OUT.stat().st_size // 1000} kB')


if __name__ == '__main__':
    main()
