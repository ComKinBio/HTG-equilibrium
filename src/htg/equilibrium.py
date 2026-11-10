"""Constrained Gibbs-energy minimization at one temperature and pressure.

    min_n  G(n) = sum over phases G_phase(n_phase)
    s.t.   elemental balances, aqueous electroneutrality,
           optional carbon-allocation constraints, n >= lower bounds

Every (phase, species) pair is an independent variable in mol per mol of
dry-feed carbon. A feasibility LP gives the starting point; the selected
optimizer (default IPOPT) minimizes G. A converged result is a local minimum
of the selected-phase problem; it is not certified to be the global one.
"""
from dataclasses import dataclass

import numpy as np
from scipy.linalg import qr
from scipy.optimize import linprog, minimize

from . import optimizer as opt
from . import standard_states as ss
from .constants import ATOMIC_WEIGHTS, PHASES
from .phases import EOSFailure, ElectrolytePhase, FluidPhase, SolidPhase, build_pcsaft_model
from .results import Result

EOS_FAILURE_TEXT = ('Invalid PC-SAFT fugacities', 'Invalid ePC-SAFT fugacity coefficients',
                    'Nonfinite phase Gibbs energy')


class SolveError(RuntimeError):
    """The calculation failed; the message starts with the stage that failed."""


def species_to_wt_percent(registry, species_moles, atomic_weights):
    """Elemental wt% of a dry feed given as relative molar amounts of registered species."""
    entries = {**registry.all_molecular, **registry.solid, **registry.ion}
    atoms = {}
    for name, amount in species_moles.items():
        if name not in entries or not np.isfinite(amount) or amount < 0:
            raise ValueError(f'Invalid dry-feed species or amount: {name}={amount}')
        for element, count in entries[name]['elements'].items():
            atoms[element] = atoms.get(element, 0.0) + amount * count
    masses = {e: count * atomic_weights[e] for e, count in atoms.items()}
    total = sum(masses.values())
    if total <= 0 or masses.get('C', 0) <= 0:
        raise ValueError('Dry feed must contain positive mass and carbon')
    return {e: mass / total * 100 for e, mass in masses.items()}


def apply_phase_rules(case, registry):
    """Effective phases, bounds, start and initialization flags at the case temperature.

    The aqueous cut-off is a modelling rule (not a computed phase boundary): at
    T >= aqueous_cutoff_C, or with aqueous_enabled=False, the aqueous phase and
    its bounds are removed. Feed atoms are kept and redistributed.
    """
    phases = {phase: list(case.phases.get(phase, [])) for phase in PHASES}
    disabled = case.aqueous_enabled is False
    cutoff = case.aqueous_cutoff_C
    if cutoff is not None:
        cutoff = float(cutoff)
        if not np.isfinite(cutoff) or cutoff <= -273.15:
            raise ValueError('aqueous_cutoff_C must be finite and above absolute zero')
        disabled = disabled or float(case.temperature_C) >= cutoff
    liquid = case.initialization in ('auto', 'liquid')
    gas_water = case.initialization == 'vapor' or (case.initialization == 'auto' and disabled)
    minimum = dict(case.minimum_amounts)
    removed = list(phases['aqueous']) if disabled else []
    if disabled:
        phases['aqueous'] = []
        minimum = {pair: v for pair, v in minimum.items() if pair[0] != 'aqueous'}
        liquid = False
        if float(case.phase_carbon_fractions.get('aqueous', 0.0)) != 0.0:
            raise ValueError('Aqueous cut-off conflicts with nonzero aqueous phase_carbon_fractions; '
                             'revise carbon targets explicitly')
    for phase, names in phases.items():
        if phase != 'aqueous':
            forbidden = [s for s in names if s in registry.ion]
            if forbidden:
                raise ValueError(f'Ions are only allowed in aqueous: {phase}: {forbidden}')
    report = {'aqueous_disabled': disabled, 'aqueous_cutoff_C': cutoff,
              'removed_aqueous_species': removed, 'initial_composition_reset': False}
    initial = case.initial_amounts
    if initial is not None:
        if any(not np.isfinite(v) or v < 0 for v in initial.values()):
            raise ValueError('Invalid initial amount.')
        selected = {(p, s) for p, names in phases.items() for s in names}
        lost = [pair for pair, v in initial.items() if pair not in selected and v > 0]
        if lost:
            initial = None
            report['initial_composition_reset'] = True
            report['initial_reset_reason'] = f'Populated species removed: {lost}'
        else:
            initial = {pair: v for pair, v in initial.items() if pair in selected}
    return phases, minimum, initial, liquid, gas_water, report


class Stage:
    """Name of the calculation stage, reported when something fails."""
    def __init__(self):
        self.name = 'case definition'


@dataclass
class Setup:
    """Everything needed to evaluate G and run the optimizer for one case."""
    case: object
    pairs: list
    phase_species: dict
    indices: dict
    phase_models: dict
    phase_G: object
    total_G: object
    initial_amounts: np.ndarray
    initial_G: float
    x_start: np.ndarray
    free: np.ndarray
    bounds: list
    A_opt: np.ndarray
    b_opt: np.ndarray
    A_full: np.ndarray
    b_full: np.ndarray
    A: np.ndarray
    b: np.ndarray
    charge: np.ndarray
    gas_carbon: np.ndarray
    feed_carbon: float
    gas_min_row: np.ndarray
    gas_min_target: float
    gas_min_fraction: float | None
    energy_scale: float
    elements: list
    entries: list
    C_mol_per_g: float
    T: float
    atomic_weights: dict
    phase_rules: dict

    def expand(self, x):
        amounts = np.zeros(len(self.pairs))
        amounts[self.free] = x
        return amounts


def solve(model, case):
    """Solve ``case`` (a single temperature) with ``model``; returns a Result."""
    stage = Stage()
    try:
        s = setup(model, case, stage)
        stage.name = 'optimization / final validation'
        return _optimize(model, s)
    except Exception as exc:
        if isinstance(exc, SolveError):
            raise
        raise SolveError(f'{stage.name}: {type(exc).__name__}: {exc}') from exc


def setup(model, case, stage=None, build_phases=True):
    """Constraints, feasible start and phase models for ``case`` (stops before optimizing).

    With build_phases=False it stops after the feasibility LP (no equation of state)
    and returns None.
    """
    stage = stage or Stage()
    case.validate()
    if (case.points is not None or isinstance(case.temperature_C, (list, tuple))
            or isinstance(case.pressure_bar, (list, tuple))):
        raise ValueError('A single temperature and pressure is needed; use Model.sweep for a grid, '
                         'or case.at(T, P).')
    registry, nasa_db = model.registry, model.databases.nasa
    stage.name = 'feed and initial composition'
    AW = {**ATOMIC_WEIGHTS, **case.atomic_weights}
    phase_species, minimum_amounts, initial_input, liquid_init, gas_water_init, allocation = \
        apply_phase_rules(case, registry)
    TEMPERATURE_C = float(case.temperature_C)
    PRESSURE_BAR = float(case.pressure_bar)
    FEED = case.feed.dry_element_wt_percent(
        lambda moles: species_to_wt_percent(registry, moles, AW))
    WATER = case.feed.total_water_g_per_g_dry()
    SALTS = case.feed.salts
    GAS_C = case.gas_carbon_fraction
    GAS_C_MIN = case.gas_carbon_min_fraction
    if GAS_C_MIN is not None:
        GAS_C_MIN = float(GAS_C_MIN)
        if not np.isfinite(GAS_C_MIN) or not 0 <= GAS_C_MIN <= 1:
            raise ValueError('Minimum gas-carbon fraction must be in [0, 1]')
        if GAS_C is not None:
            raise ValueError('Use either gas_carbon_fraction (equality) or '
                             'gas_carbon_min_fraction (lower bound), not both')
    SOLID_C = case.solid_carbon_fraction
    WATER_FLOOR = float(case.aqueous_water_floor)
    AQUEOUS_MODEL = case.aqueous_model
    if AQUEOUS_MODEL == 'molecular' and any(s in registry.ion for s in phase_species['aqueous']):
        raise ValueError('Ions require the electrolyte aqueous model')
    INITIAL_GAS_WATER_FRACTION = float(case.initial_gas_water_fraction)
    if not np.isfinite(INITIAL_GAS_WATER_FRACTION) or INITIAL_GAS_WATER_FRACTION < 0:
        raise ValueError('Invalid initial gas-water fraction')
    if not np.isfinite(TEMPERATURE_C) or TEMPERATURE_C <= -273.15:
        raise ValueError('Invalid temperature.')
    if not np.isfinite(PRESSURE_BAR) or PRESSURE_BAR <= 0:
        raise ValueError('Pressure must be finite and positive.')
    if not np.isfinite(WATER) or WATER < 0:
        raise ValueError('Added water must be finite and nonnegative.')
    if any(not np.isfinite(v) or v < 0 for v in FEED.values()):
        raise ValueError('Invalid feed elemental composition.')
    if FEED.get('C', 0.0) <= 0:
        raise ValueError('This carbon-normalized calculation requires positive dry-feed carbon.')
    if sum(FEED.values()) > 100.0 + 1e-08:
        raise ValueError('Feed elemental percentages exceed 100%.')
    if phase_species['aqueous'] and AQUEOUS_MODEL == 'electrolyte' and (
            not np.isfinite(WATER_FLOOR) or WATER_FLOOR <= 1e-09):
        raise ValueError('aqueous_water_floor must exceed the aqueous model\'s 1e-9 mol '
                         'solvent cutoff.')
    for phase, names in phase_species.items():
        if len(names) != len(set(names)):
            raise ValueError(f'Duplicate species in {phase}.')

    pairs = [(phase, s) for phase, names in phase_species.items() for s in names]
    if not pairs:
        raise ValueError('Select at least one species.')
    entries = [registry.entry(phase, s) for phase, s in pairs]
    indices = {phase: np.array([j for j, (p, _) in enumerate(pairs) if p == phase], dtype=int)
               for phase in PHASES}
    elements = sorted(set(FEED) | {'H', 'O'}
                      | {e for salt in SALTS.values() for e in salt.elements}
                      | {e for entry in entries for e in entry['elements']})
    missing_aw = [e for e in elements if e not in AW]
    if missing_aw:
        raise ValueError('Missing atomic weights: ' + ', '.join(missing_aw))
    A = np.array([[entry['elements'].get(e, 0.0) for entry in entries] for e in elements],
                 dtype=float)
    C_mol_per_g = FEED['C'] / 100.0 / AW['C']
    b = np.array([FEED.get(e, 0.0) / 100.0 / AW[e] / C_mol_per_g for e in elements], dtype=float)
    water_MW = 2.0 * AW['H'] + AW['O']
    water_mol_per_C = WATER / water_MW / C_mol_per_g
    b[elements.index('H')] += 2.0 * water_mol_per_C
    b[elements.index('O')] += water_mol_per_C
    for salt_name, salt in SALTS.items():
        amount = float(salt.mol_per_g_dry_feed)
        if not np.isfinite(amount) or amount < 0:
            raise ValueError(f'Invalid amount for {salt_name}.')
        for element, count in salt.elements.items():
            if not np.isfinite(count) or count < 0:
                raise ValueError(f'Invalid elemental count for {salt_name}.')
            b[elements.index(element)] += amount * count / C_mol_per_g
    carbon = A[elements.index('C')]
    feed_carbon = b[elements.index('C')]
    gas_carbon = carbon * np.array([phase == 'gas' for phase, _ in pairs], dtype=float)
    solid_carbon = np.zeros(len(pairs))
    if SOLID_C is not None:
        if case.solid_carbon_target == 'char':
            if ('solid', 'carbon') not in pairs:
                raise ValueError("Select solid 'carbon' to constrain char carbon.")
            j = pairs.index(('solid', 'carbon'))
            solid_carbon[j] = carbon[j]
        else:
            solid_carbon = carbon * np.array([phase == 'solid' for phase, _ in pairs], dtype=float)
    charge = np.array([registry.ion[s]['charge'] if phase == 'aqueous' and s in registry.ion
                       else 0.0 for phase, s in pairs], dtype=float)
    missing = [e for e, bi, row in zip(elements, b, A) if bi > 0 and not np.any(row > 0)]
    if missing:
        raise ValueError('No allowed product carrier for feed elements: ' + ', '.join(missing)
                         + '. Add neutral gas species or solids that contain them.')
    rows, targets, carbon_constraints = list(A) + [charge], list(b) + [0.0], []
    if GAS_C is not None:
        fraction = float(GAS_C)
        if not np.isfinite(fraction) or not 0 <= fraction <= 1:
            raise ValueError('Gas-carbon fraction must be in [0, 1].')
        rows.append(gas_carbon)
        targets.append(fraction * feed_carbon)
        carbon_constraints.append((gas_carbon, fraction))
    if SOLID_C is not None:
        fraction = float(SOLID_C)
        if not np.isfinite(fraction) or not 0 <= fraction <= 1:
            raise ValueError('Solid-carbon fraction must be in [0, 1].')
        if GAS_C is not None and float(GAS_C) + fraction > 1.0 + 1e-10:
            raise ValueError('Gas-carbon plus solid-carbon fractions exceed one.')
        rows.append(solid_carbon)
        targets.append(fraction * feed_carbon)
        carbon_constraints.append((solid_carbon, fraction))
    phase_fractions = case.phase_carbon_fractions
    if set(phase_fractions) - set(PHASES):
        raise ValueError('Unknown phase in phase_carbon_fractions')
    fractions = {phase: float(v) for phase, v in phase_fractions.items()}
    if any(not np.isfinite(v) or not 0 <= v <= 1 for v in fractions.values()):
        raise ValueError('Phase-carbon fractions must be finite and in [0, 1]')
    if sum(fractions.values()) > 1 + 1e-10:
        raise ValueError('Phase-carbon fractions exceed one')
    if 'gas' in fractions and GAS_C is not None and not np.isclose(
            fractions['gas'], GAS_C, atol=1e-10, rtol=0):
        raise ValueError('Conflicting gas-carbon constraints')
    for phase, fraction in fractions.items():
        row = carbon * np.array([p == phase for p, _ in pairs], dtype=float)
        rows.append(row)
        targets.append(fraction * feed_carbon)
        carbon_constraints.append((row, fraction))
    A_full = np.asarray(rows, dtype=float)
    b_full = np.asarray(targets, dtype=float)
    lower = np.zeros(len(pairs))
    for pair, value in minimum_amounts.items():
        if pair not in pairs:
            raise ValueError(f'Lower-bound species not selected: {pair}')
        if not np.isfinite(value) or value < 0:
            raise ValueError(f'Invalid lower bound for {pair}.')
        lower[pairs.index(pair)] = value
    if phase_species['aqueous']:
        if 'water' not in phase_species['aqueous']:
            raise ValueError('The aqueous phase requires solvent water.')
        w = pairs.index(('aqueous', 'water'))
        lower[w] = max(lower[w], WATER_FLOOR)
    fixed_zero = np.zeros(len(pairs), dtype=bool)
    # Species containing an element absent from the feed must be zero.
    for row, inventory in zip(A, b):
        if inventory == 0.0:
            fixed_zero |= row > 0
    for row, fraction in carbon_constraints:
        if fraction == 0.0:
            fixed_zero |= row > 0
        elif fraction == 1.0:
            fixed_zero |= (carbon > 0) & (row == 0)
    if np.any(lower[fixed_zero] > 0):
        raise ValueError('A minimum amount conflicts with a carbon constraint.')
    free = np.flatnonzero(~fixed_zero)
    if not len(free):
        raise ValueError('No free species remain.')

    def expand(x):
        amounts = np.zeros(len(pairs))
        amounts[free] = x
        return amounts

    reduced = A_full[:, free]
    _, R_qr, pivots = qr(reduced.T, pivoting=True, mode='economic')
    diagonal = np.abs(np.diag(R_qr))
    rank = int(np.sum(diagonal > 1e-10 * max(1.0, float(diagonal.max()) if diagonal.size else 0.0)))
    selected_rows = np.sort(pivots[:rank])
    A_opt = reduced[selected_rows]
    b_opt = b_full[selected_rows]
    bounds = [(float(lower[j]), None) for j in free]
    start_objective = np.zeros(len(free))
    if gas_water_init:
        for k, j in enumerate(free):
            if pairs[j] == ('gas', 'water'):
                start_objective[k] = -1.0
    start_bounds = list(bounds)
    if liquid_init and ('aqueous', 'water') in pairs:
        for k, j in enumerate(free):
            if pairs[j] == ('aqueous', 'water'):
                start_objective[k] = -1.0
            elif pairs[j] == ('gas', 'water'):
                lower_bound, _ = start_bounds[k]
                start_bounds[k] = (lower_bound, max(lower_bound, INITIAL_GAS_WATER_FRACTION
                                                    * water_mol_per_C))
    # The gas-carbon minimum is an inequality in the LP, optimizer and final check.
    gas_min_row = gas_carbon[free]
    gas_min_target = (GAS_C_MIN or 0.0) * feed_carbon
    gas_min_lp = ({'A_ub': -gas_min_row[None, :], 'b_ub': np.array([-gas_min_target])}
                  if GAS_C_MIN is not None else {})
    feasible = linprog(start_objective, A_eq=A_opt, b_eq=b_opt, bounds=start_bounds,
                       method='highs', **gas_min_lp)
    if not feasible.success:
        raise ValueError('Selected feed, species and constraints are infeasible: '
                         + feasible.message)
    initial_amounts = expand(feasible.x)
    if initial_input is not None:
        initial_amounts = np.array([initial_input.get(pair, 0.0) for pair in pairs], dtype=float)
        if set(initial_input) - set(pairs):
            raise ValueError('Initial composition contains an unselected phase/species')
        if (np.any(~np.isfinite(initial_amounts)) or np.any(initial_amounts < lower)
                or np.any(initial_amounts[fixed_zero] != 0)):
            raise ValueError('Initial composition violates species bounds')
        feasible.x = initial_amounts[free].copy()
    if np.max(np.abs(A_full @ initial_amounts - b_full)) > 1e-08:
        raise ValueError('The selected species cannot satisfy all original balances and '
                         'carbon constraints.')
    if GAS_C_MIN is not None and gas_carbon @ initial_amounts < gas_min_target - 1e-8:
        raise ValueError('Initial composition violates minimum gas-carbon fraction')

    if not build_phases:
        return None
    stage.name = 'thermochemistry / aqueous validity'
    T = TEMPERATURE_C + 273.15
    P = PRESSURE_BAR * 100000.0
    electrolyte = bool(phase_species['aqueous']) and AQUEOUS_MODEL == 'electrolyte'
    to_resolve = [s for phase, names in phase_species.items() if phase != 'solid'
                  for s in names if s not in registry.ion]
    if electrolyte:
        to_resolve += [s for s in registry.aqueous_reference_neutrals]
    resolved = registry.resolve_parameters(model.databases, T, list(dict.fromkeys(to_resolve)), AW)
    phase_models = {}
    if electrolyte:
        phase_models['aqueous'] = ElectrolytePhase(model, phase_species['aqueous'], T, P,
                                                   resolved, AW)
    for phase, names in phase_species.items():
        if not names or phase in phase_models:
            continue
        if phase == 'solid':
            phase_models[phase] = SolidPhase(names, ss.solid_g0(registry, nasa_db, names, T, P))
        else:
            root = 'vapour' if phase == 'gas' else 'liquid'
            phase_models[phase] = FluidPhase(
                phase, names, ss.nasa_g0(registry, nasa_db, names, T),
                build_pcsaft_model(names, resolved, registry.binary_interactions), T, P, root)
    phase_models = {phase: phase_models[phase] for phase in PHASES if phase in phase_models}

    def phase_G(amounts):
        terms = {}
        for phase in PHASES:
            n = amounts[indices[phase]]
            if not phase_species[phase] or n.sum() <= 0.0:
                terms[phase] = 0.0
            else:
                terms[phase] = float(phase_models[phase].gibbs(n))
        if not np.all(np.isfinite(list(terms.values()))):
            raise ValueError('Nonfinite phase Gibbs energy.')
        return terms

    def total_G(amounts):
        if np.min(amounts) < -1e-10:
            raise ValueError('Negative species amount.')
        return float(sum(phase_G(np.maximum(amounts, 0.0)).values()))

    stage.name = 'initial Gibbs evaluation'
    initial_G = total_G(initial_amounts)
    energy_scale = float(case.energy_scale if case.energy_scale is not None
                         else max(10000.0, abs(initial_G)))
    if not np.isfinite(energy_scale) or energy_scale <= 0:
        raise ValueError('energy_scale must be finite and positive')
    return Setup(
        case=case, pairs=pairs, phase_species=phase_species, indices=indices,
        phase_models=phase_models, phase_G=phase_G, total_G=total_G,
        initial_amounts=initial_amounts, initial_G=initial_G, x_start=feasible.x, free=free,
        bounds=bounds, A_opt=A_opt, b_opt=b_opt, A_full=A_full, b_full=b_full, A=A, b=b,
        charge=charge, gas_carbon=gas_carbon, feed_carbon=feed_carbon, gas_min_row=gas_min_row,
        gas_min_target=gas_min_target, gas_min_fraction=GAS_C_MIN, energy_scale=energy_scale,
        elements=elements, entries=entries, C_mol_per_g=C_mol_per_g, T=T, atomic_weights=AW,
        phase_rules=allocation)


def _optimize(model, s):
    case, registry = s.case, model.registry
    pairs, expand, total_G, phase_G = s.pairs, s.expand, s.total_G, s.phase_G
    A_opt, b_opt, bounds, initial_G = s.A_opt, s.b_opt, s.bounds, s.initial_G
    gas_min_row, gas_min_target, GAS_C_MIN = s.gas_min_row, s.gas_min_target, s.gas_min_fraction
    gas_carbon, charge, energy_scale = s.gas_carbon, s.charge, s.energy_scale
    x_start = s.x_start

    eos_errors = []

    constraints = [{'type': 'eq', 'fun': lambda x: A_opt @ x - b_opt, 'jac': lambda x: A_opt}]
    if GAS_C_MIN is not None:
        constraints.append({'type': 'ineq', 'fun': lambda x: gas_min_row @ x - gas_min_target,
                            'jac': lambda x: gas_min_row})
    if case.solver == 'slsqp-fd':
        # The original method: SLSQP with finite-difference gradients.
        def objective(x):
            try:
                return (total_G(expand(x)) - initial_G) / energy_scale
            except ValueError as exc:
                message = str(exc)
                if not any(text in message for text in EOS_FAILURE_TEXT):
                    raise
                if len(eos_errors) < 10:
                    eos_errors.append(message)
                return float(1000000.0 + np.sum((x - x_start) ** 2))
        optimized = minimize(objective, x_start, method='SLSQP', bounds=bounds,
                             constraints=constraints,
                             options={'maxiter': int(case.max_iterations),
                                      'ftol': float(case.fd_ftol), 'eps': float(case.fd_eps),
                                      'disp': False})
        diagnostics = {'method': 'slsqp-fd'}
    else:
        use_min = GAS_C_MIN is not None
        problem = opt.GibbsProblem(
            phases=s.phase_models, pairs=pairs, free=s.free, phase_indices=s.indices,
            ions=set(registry.ion), T=s.T, A_eq=A_opt, b_eq=b_opt,
            ineq_row=gas_min_row if use_min else None,
            ineq_target=gas_min_target if use_min else None,
            lower=[low for low, _ in bounds], x_start=x_start, G0=initial_G,
            scale=energy_scale, model_total_G=total_G)
        options = {'method': case.solver, 'max_iterations': int(case.max_iterations),
                   'multistart': int(case.multistart), 'seed': int(case.seed),
                   **case.solver_options}
        if case.stability_test is not None:
            options['stability_test'] = bool(case.stability_test)
        optimized, diagnostics = opt.optimize(problem, options)
        eos_errors.extend(problem.eos_errors[:max(0, 10 - len(eos_errors))])
    amounts = expand(optimized.x)
    residual = float(np.max(np.abs(s.A_full @ amounts - s.b_full)))
    gas_min_violation = (max(0.0, float(gas_min_target - gas_carbon @ amounts))
                         if GAS_C_MIN is not None else 0.0)
    if gas_min_violation > 1e-7:
        raise RuntimeError(f'Minimum gas-carbon allocation violated by {gas_min_violation:.3e} '
                           'mol C/mol dry-feed C')
    if (not optimized.success or np.any(~np.isfinite(amounts)) or np.min(amounts) < -1e-09
            or residual > 1e-07):
        raise RuntimeError(f'{optimized.message}; constraint residual={residual:.3e}; '
                           f'recorded EOS errors={eos_errors}')
    amounts = np.maximum(amounts, 0.0)
    G_terms = phase_G(amounts)
    final_G = float(sum(G_terms.values()))
    if final_G > initial_G + 0.0001:
        raise RuntimeError('Final Gibbs energy exceeds the feasible start.')
    return Result(
        temperature_C=float(case.temperature_C), pressure_bar=float(case.pressure_bar),
        amounts=dict(zip(pairs, map(float, amounts))),
        initial_amounts=dict(zip(pairs, map(float, s.initial_amounts))),
        gibbs_J_per_mol_C=final_G, phase_gibbs_J_per_mol_C=G_terms,
        balance_residual=residual, charge_residual=float(charge @ amounts),
        iterations=int(optimized.nit), eos_trial_errors=list(eos_errors),
        feed_C_mol_per_g_dry=s.C_mol_per_g, elements=s.elements, feed_atoms=s.b.copy(),
        product_atoms=s.A @ amounts,
        gas_carbon_fraction=float(gas_carbon @ amounts / s.feed_carbon),
        gas_carbon_min_violation=gas_min_violation, phase_rules=s.phase_rules,
        optimizer=diagnostics, atomic_weights=s.atomic_weights,
        species_data={pair: {'formula': entry.get('formula', pair[1]),
                             'elements': dict(entry['elements']),
                             'charge': entry.get('charge', 0) if pair[1] in registry.ion else 0}
                      for pair, entry in zip(pairs, s.entries)})
