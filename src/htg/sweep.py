"""Temperature sweeps and preflight checks."""
import pandas as pd

from . import standard_states as ss
from .equilibrium import apply_phase_rules, setup
from .results import Result, SweepResult
from .runner import WorkerPool, solve_chain


def run_sweep(model, case, progress=print, workers=1):
    """Solve ``case`` at every point of its temperature-pressure grid.

    With continuation each chain (rising temperature at one pressure, or the
    order of explicit points) starts each point from the previous successful
    solution when the selected (phase, species) pairs are unchanged; otherwise
    from a fresh feasible LP. A failed point is recorded and the sweep continues.
    A heating path can follow a local branch; smooth profiles do not prove
    global equilibrium.

    workers > 1 solves independent chains in parallel worker processes (each
    starts its own Julia, which takes about a minute).
    """
    progress = progress or (lambda message: None)
    chains = case.chains()
    if workers > 1 and len(chains) > 1:
        return _parallel(model, case, progress, min(workers, len(chains)))
    outcomes = {}

    def on_point(T, P, outcome, start):
        outcomes[(T, P)] = outcome
        text = 'solved' if isinstance(outcome, Result) else f'FAILED: {outcome}'
        progress(f'{T:g} °C, {P:g} bar: {text} (start: {start})')

    for chain in chains:
        solve_chain(model, case, chain, on_point)
    return SweepResult(case, [(T, P, outcomes[(T, P)]) for T, P in case.grid()])


def _parallel(model, case, progress, workers):
    progress(f'Solving {len(case.grid())} points in {workers} worker processes '
             '(each starts Julia; the first results take about a minute)')

    def on_point(run, T, P, outcome, start):
        text = 'solved' if isinstance(outcome, Result) else f'FAILED: {outcome}'
        progress(f'{T:g} °C, {P:g} bar: {text} (start: {start})')

    with WorkerPool(workers, model.load_args, warm_up=False, on_point=on_point) as pool:
        run = pool.submit(case)
        run.wait()
    return run.sweep_result()


def preflight(model, case):
    """Check every point without the equation of state.

    Verifies species names, standard-state and SAFT data coverage at each
    temperature, aqueous validity limits, and feasibility of the full set of
    balances and constraints. Passing does not guarantee that the solver converges.
    """
    registry, nasa_db = model.registry, model.databases.nasa
    rows = []
    for T, P_bar in case.grid():
        point = case.at(T, P_bar)
        row = {'T_C': T, 'P_bar': P_bar, 'status': 'ready', 'aqueous_removed': None,
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
