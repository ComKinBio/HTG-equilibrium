"""Temperature sweeps and preflight checks."""
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor

import pandas as pd

from . import standard_states as ss
from .equilibrium import SolveError, apply_phase_rules, setup
from .results import Result, SweepResult


def run_sweep(model, case, progress=print, workers=1):
    """Solve ``case`` at every temperature of its grid (sorted, with refinement points).

    continuation=True: each point starts from the last successful solution when
    the selected (phase, species) pairs are unchanged; otherwise from a fresh
    feasible LP. A failed point is recorded and the sweep continues. A heating
    path can follow a local branch; smooth profiles do not prove global equilibrium.

    continuation=False and workers > 1 solves points in parallel processes
    (each process starts its own Julia, so this pays off only for long sweeps).
    """
    progress = progress or (lambda message: None)
    grid = case.temperatures()
    if not case.continuation and workers > 1 and len(grid) > 1:
        return _parallel(model, case, grid, progress, workers)
    points, previous = [], None
    for T in grid:
        point = case.at(T)
        start = 'feasible LP'
        if case.continuation and previous is not None:
            phases = apply_phase_rules(point, model.registry)[0]
            signature = {(p, s) for p, names in phases.items() for s in names}
            if set(previous.amounts) == signature:
                point.initial_amounts = dict(previous.amounts)
                start = f'solution at {previous.temperature_C:g} °C'
            else:
                point.initial_amounts = None
                start = 'feasible LP (phase selection changed)'
        try:
            result = model.solve(point)
            points.append((T, result))
            previous = result
            progress(f'{T:g} °C: solved (start: {start})')
        except SolveError as exc:
            points.append((T, str(exc)))
            progress(f'{T:g} °C: FAILED (start: {start}): {exc}')
    return SweepResult(case, points)


_worker_model = None


def _worker_init(data_dir, nasa_names):
    global _worker_model
    from .model import Model
    _worker_model = Model.load(data_dir, nasa_names=nasa_names)


def _worker_solve(point):
    try:
        return _worker_model.solve(point)
    except SolveError as exc:
        return str(exc)


def _parallel(model, case, grid, progress, workers):
    workers = min(len(grid), workers or max(1, (os.cpu_count() or 2) // 2))
    progress(f'Solving {len(grid)} temperatures in {workers} processes '
             '(each starts Julia; the first results take a while)')
    nasa_names = model.registry.config.get('nasa_names', {})
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn'),
                             initializer=_worker_init, initargs=(str(model.data_dir), nasa_names)) as pool:
        outcomes = list(pool.map(_worker_solve, [case.at(T) for T in grid]))
    for T, outcome in zip(grid, outcomes):
        progress(f'{T:g} °C: ' + ('solved' if isinstance(outcome, Result) else f'FAILED: {outcome}'))
    return SweepResult(case, list(zip(grid, outcomes)))


def preflight(model, case):
    """Check every temperature without the equation of state.

    Verifies species names, standard-state and SAFT data coverage at each
    temperature, aqueous validity limits, and feasibility of the full set of
    balances and constraints. Passing does not guarantee that the solver converges.
    """
    registry, nasa_db = model.registry, model.databases.nasa
    rows = []
    for T in case.temperatures():
        point = case.at(T)
        row = {'T_C': T, 'P_bar': point.pressure_bar, 'status': 'ready', 'aqueous_removed': None,
               'error': ''}
        try:
            phases, *_, report = apply_phase_rules(point, registry)
            row['aqueous_removed'] = report['aqueous_disabled']
            setup(model, point, build_phases=False)
            T_K, P = T + 273.15, point.pressure_bar * 1e5
            for phase, names in phases.items():
                for name in names:
                    if phase == 'solid':
                        ss.solid_g0(registry, nasa_db, [name], T_K, P)
                    elif name not in registry.ion:
                        ss.nasa_g0(registry, nasa_db, [name], T_K)
                        registry.resolve_parameters(model.databases, T_K, [name])
            if phases['aqueous'] and point.aqueous_model == 'electrolyte':
                water = ss.supcrt_water(T_K, P)
                if not (water['rho'] >= 350 and water['epsilon'] > 1):
                    raise ValueError(f"Aqueous electrolyte model invalid here: water density "
                                     f"{water['rho']:.0f} kg/m3 (needs >= 350). Set aqueous "
                                     "cutoff_C below this temperature.")
        except Exception as exc:
            row.update(status='invalid', error=f'{type(exc).__name__}: {exc}')
        rows.append(row)
    return pd.DataFrame(rows)
