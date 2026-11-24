"""Dump the Julia-side inputs of every solved reference point (M1 spike).

For each solved case in tests/reference/data/solves.json (run with the app's
default solver, IPOPT) this records exactly what the Python app hands to
Julia for one point:

- the Clapeyron model inputs (PCSAFT components + userlocations; for the
  electrolyte phase also the ePCSAFT/hsdDH/ConstRSP arguments),
- the phase data computed in Python (standard-state mu0, HKF ion g0, T, P,
  root, probe weights, positions),
- the IPOPT problem (A, b, inequalities, lower bounds, starts, G0, scale,
  options) exactly as optimizer._ipopt builds it,
- check values computed by Python through juliacall (fugacity coefficients
  and objective/gradient/Hessian at the start, ln phi at infinite dilution),
- the Python app's own raw IPOPT result from every start and the full solve
  result (G, amounts), for comparison with the pure-Julia replay.

Run from the repository root:
    uv run python julia/spikes/m1_threads/dump_problems.py
Writes julia/spikes/m1_threads/problems.json. Does not change the package.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / 'tests' / 'reference'))

import htg  # noqa: E402
import pyclapeyron  # noqa: E402
from htg import optimizer as opt  # noqa: E402
from htg.equilibrium import setup  # noqa: E402
from htg.phases import ElectrolytePhase, FluidPhase, SolidPhase  # noqa: E402

from convert import case_from_legacy, load  # noqa: E402

# Record every PCSAFT construction (the only Clapeyron constructor called with
# userlocations); the phase objects keep the returned model, so identity works.
_RECORDS = {}
_ORIGINAL_PCSAFT = pyclapeyron.PCSAFT


def _recording_pcsaft(components, userlocations=None, **kwargs):
    model = _ORIGINAL_PCSAFT(components, userlocations=userlocations, **kwargs)
    _RECORDS[id(model)] = (model, list(components), userlocations, kwargs)
    return model


pyclapeyron.PCSAFT = _recording_pcsaft


def _floats(a):
    return [float(v) for v in np.asarray(a, dtype=float).reshape(-1)]


def _userlocations(model):
    _, components, ul, kwargs = _RECORDS[id(model)]
    out = {'components': components, 'verbose': kwargs.get('verbose', None)}
    for key in ('Mw', 'segment', 'sigma', 'epsilon'):
        out[key] = _floats(ul[key])
    for key in ('n_H', 'n_e'):
        if key in ul:
            out[key] = [int(v) for v in ul[key]]
    for key in ('epsilon_assoc', 'bondvol'):
        if key in ul:
            out[key] = [[list(a), list(b), float(v)] for (a, b), v in ul[key].items()]
    if 'k' in ul:
        out['k'] = [_floats(row) for row in np.asarray(ul['k'])]
    return out


def _phase_record(name, phase, x_amounts):
    if isinstance(phase, SolidPhase):
        return {'type': 'solid', 'species': phase.species, 'mu0': _floats(phase.mu0)}
    if isinstance(phase, FluidPhase):
        rec = {'type': 'fluid', 'name': phase.name, 'species': phase.species, 'mu0': _floats(phase.mu0),
               'T': phase.T, 'P': phase.P, 'root': phase.root, 'pcsaft': _userlocations(phase.eos_model)}
        positive = np.maximum(x_amounts, 1e-14)
        rec['check_n'] = _floats(positive)
        rec['check_phi'] = _floats(phase._phis(positive))
        return rec
    if isinstance(phase, ElectrolytePhase):
        rec = {'type': 'electrolyte', 'species': phase.species, 'T': phase.T, 'P': phase.P,
               'all_species': phase.all_species, 'neutrals': phase.neutrals, 'ions': phase.ions,
               'positions': [int(p) for p in phase.positions], 'z': _floats(phase.z),
               'probe_weights': [float(phase.probe_weights.get(s, 1.)) for s in phase.ions],
               'components': list(phase.components), 'charges': [int(c) for c in phase.charges],
               'epsilon_r': float(phase.water['epsilon']),
               'short_range': _userlocations(phase.short_range),
               'mu0_ion': _floats(phase.mu0_ion), 'mu0_neutral': _floats(phase.mu0_neutral),
               'lnphi_infinity': _floats(phase.lnphi_infinity)}
        full = phase.expand(x_amounts)
        if full[0] > 1e-9:
            rec['check_n'] = _floats(full)
            rec['check_phi'] = _floats(phase.phi(full))
        return rec
    raise TypeError(type(phase))


def _ipopt_options(options):
    """The option dict optimizer._ipopt passes to Ipopt (same order and kinds)."""
    opts = {'tol': 1e-10, 'acceptable_tol': 1e-7, 'max_iter': int(options.get('max_iterations', 3000)),
            'print_level': int(options.get('print_level', 0)),
            'mu_strategy': 'monotone', 'mu_init': 1e-9,
            'bound_relax_factor': 0.0, 'honor_original_bounds': 'yes',
            'bound_push': 1e-6, 'bound_frac': 1e-6,
            'jac_c_constant': 'yes', 'jac_d_constant': 'yes', 'sb': 'yes'}
    opts.update(options.get('ipopt_options', {}))
    return [[k, v, 'str' if isinstance(v, str) else 'int' if isinstance(v, int) else 'num']
            for k, v in opts.items()]


def dump_case(models, entry):
    case = case_from_legacy(entry)
    if case.solver != 'ipopt':
        case = case.replace(solver='ipopt')
    model = models[entry['model']]
    s = setup(model, case)
    registry = model.registry
    problem = opt.GibbsProblem(
        phases=s.phase_models, pairs=s.pairs, free=s.free, phase_indices=s.indices,
        ions=set(registry.ion), T=s.T, A_eq=s.A_opt, b_eq=s.b_opt,
        ineq_rows=s.ineq_rows if len(s.ineq_rows) else None,
        ineq_targets=s.ineq_targets if len(s.ineq_rows) else None,
        lower=[low for low, _ in s.bounds], x_start=s.x_start, G0=s.initial_G,
        scale=s.energy_scale, model_total_G=s.total_G)
    options = {'method': case.solver, 'max_iterations': int(case.max_iterations),
               'starts': int(case.starts), 'seed': int(case.seed), **case.solver_options}
    amounts0 = problem.expand(problem.x_start)
    phases = {name: _phase_record(name, ph, amounts0[problem.phase_indices[name]])
              for name, ph in problem.phases.items()}
    f0, g0 = problem.fun_and_grad(problem.x_start, penalty=False)
    hblocks = {ph: {'pos': [int(p) for p in pos], 'H': [_floats(r) for r in H]}
               for ph, (pos, H) in problem.hessian_blocks(problem.x_start).items()}
    starts = opt.candidate_starts(problem, int(case.starts), int(case.seed))
    python_runs = []
    for label, x0 in starts:
        problem.nfev = 0
        t0 = time.perf_counter()
        res = opt._ipopt(problem, x0, options)
        python_runs.append({'start': label, 'x0': _floats(x0), 'x': _floats(res.x), 'status': int(res.status),
                            'obj': float(res.fun), 'hessian_evals': int(res.nit), 'nfev': problem.nfev,
                            'seconds': time.perf_counter() - t0})
    t0 = time.perf_counter()
    full = model.solve(case)
    full_seconds = time.perf_counter() - t0
    ref = entry['result']
    return {
        'name': entry['name'], 'model': entry['model'], 'reference_solver': (ref.get('optimizer') or {}).get(
            'method', 'slsqp-fd'),
        'T': s.T, 'P': float(case.pressure_bar) * 1e5,
        'pairs': [list(p) for p in s.pairs], 'free': [int(j) for j in s.free],
        'phase_order': list(problem.phases), 'phase_indices': {k: [int(i) for i in problem.phase_indices[k]]
                                                               for k in problem.phases},
        'phases': phases,
        'A_eq': [_floats(r) for r in problem.A_eq], 'b_eq': _floats(problem.b_eq),
        'ineq_rows': None if problem.ineq_rows is None else [_floats(r) for r in problem.ineq_rows],
        'ineq_targets': None if problem.ineq_targets is None else _floats(problem.ineq_targets),
        'lower': _floats(problem.lower), 'x_start': _floats(problem.x_start), 'G0': problem.G0,
        'scale': problem.scale, 'ipopt_options': _ipopt_options(options),
        'check': {'f0': float(f0), 'grad0': _floats(g0), 'hessian_blocks0': hblocks},
        'python_ipopt': python_runs,
        'python_solve': {'G': full.gibbs_J_per_mol_C, 'amounts': _floats(list(full.amounts.values())),
                         'seconds': full_seconds, 'polish_used': bool(
                             full.optimizer.get('polish', {}).get('used', False))},
        'reference': {'G': ref['gibbs_J_per_mol_feed_C'],
                      'amounts': [float(v) for _, _, v in ref['amounts']]},
    }


def main():
    FROZEN_DATA = REPO / 'tests' / 'reference' / 'databases'  # as in tests/conftest.py
    models = {'base': htg.Model.load(FROZEN_DATA, nasa_names={}), 'acetone': htg.Model.load(FROZEN_DATA)}
    out = []
    for entry in load('solves')['cases']:
        if 'result' not in entry:
            continue
        t0 = time.perf_counter()
        out.append(dump_case(models, entry))
        r = out[-1]
        print(f"{r['name']:32s} n={len(r['x_start']):3d} starts={len(r['python_ipopt'])} "
              f"G={r['python_solve']['G']:.6f} ref={r['reference']['G']:.6f} "
              f"({time.perf_counter() - t0:.1f} s)", flush=True)
    with open(HERE / 'problems.json', 'w') as f:
        json.dump({'cases': out}, f)


if __name__ == '__main__':
    main()
