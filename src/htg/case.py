"""Case definition: everything that specifies one equilibrium calculation or sweep.

A case can be written in Python::

    case = htg.Case(temperature_C=450, pressure_bar=240,
                    feed=htg.Feed(elements_wt_percent={'C': 30.4, ...}, water_g_per_g_dry=4.0),
                    phases={'gas': [...], 'aqueous': [...], 'solid': ['carbon']})

or as a TOML case file (see cases/ and ``htg new``) and loaded with ``htg.load_case``.

Units: temperatures in °C, pressures in bar (absolute). Species amounts in
constraints and results are mol per mol of dry-feed carbon.
"""
import copy
import math
import tomllib
from dataclasses import MISSING, asdict, dataclass, field, fields
from pathlib import Path

import numpy as np

from .constants import PHASES

SOLVERS = ('ipopt', 'slsqp', 'trust-constr', 'slsqp-fd')
INITIALIZATIONS = ('auto', 'liquid', 'vapor', 'lp')


@dataclass
class Salt:
    """An added salt: supplies atoms only; its identity does not constrain the products."""
    mol_per_g_dry_feed: float
    elements: dict


@dataclass
class Feed:
    """Dry feed, water and salts.

    Give the dry feed either as an elemental analysis (``elements_wt_percent``)
    or as relative molar amounts of registered species (``species_moles``).

    elements_wt_percent  elemental mass percentages (ultimate analysis)
    basis                'dry': per 100 g dry feed including ash;
                         'dry_ash_free': must sum to 100, scaled by (1 - ash/100)
    ash_wt_percent       dry-basis ash, outside the reactive model (bookkeeping only)
    unmodelled_wt_percent other dry mass outside the model
                         If either of these two is given, elements + ash + unmodelled
                         must sum to 100 (within 0.05).
    species_moles        e.g. {'cellulose': 1.0}; replaces the elemental analysis
    moisture_wt_percent  feed moisture, % of as-received (wet) feed
    water_g_per_g_dry    added process water [g water / g dry feed]
    salts                {name: Salt}, in mol salt / g dry feed
    species_amounts      format-2 species amounts {name: {value, unit}} including
                         water; converted with molar masses by ``resolved``
                         (the calculation does this; see htg.feed_v2)
    """
    elements_wt_percent: dict | None = None
    basis: str = 'dry'
    ash_wt_percent: float | None = None
    unmodelled_wt_percent: float | None = None
    species_moles: dict | None = None
    moisture_wt_percent: float = 0.0
    water_g_per_g_dry: float = 0.0
    salts: dict = field(default_factory=dict)
    species_amounts: dict | None = None

    def __post_init__(self):
        self.salts = {name: s if isinstance(s, Salt) else Salt(**s) for name, s in self.salts.items()}

    def resolved(self, molar_mass):
        """This feed with species amounts converted to species_moles and a water
        ratio, using ``molar_mass(name) -> g/mol``. Other feeds are returned as is."""
        if self.species_amounts is None:
            return self
        from .feed_v2 import resolve_species
        moles, water = resolve_species(self.species_amounts, molar_mass)
        return Feed(species_moles=moles, water_g_per_g_dry=water, salts=copy.deepcopy(self.salts))

    def _unresolved(self):
        if self.species_amounts is not None:
            raise ValueError('The feed is given as species amounts; convert it with '
                             'Feed.resolved(molar_mass) first.')

    def total_water_g_per_g_dry(self):
        self._unresolved()
        moisture = _number(self.moisture_wt_percent, 'feed moisture_wt_percent', 0)
        if moisture >= 100:
            raise ValueError('Feed moisture must be below 100 wt%.')
        return _number(self.water_g_per_g_dry, 'feed water_g_per_g_dry', 0) + moisture / (100 - moisture)

    def dry_element_wt_percent(self, species_to_wt_percent):
        """Dry-basis elemental wt% used by the model.

        ``species_to_wt_percent`` converts ``species_moles`` (it needs the registry).
        """
        self._unresolved()
        ash = None if self.ash_wt_percent is None else _number(self.ash_wt_percent, 'ash_wt_percent', 0, 100)
        other = (None if self.unmodelled_wt_percent is None
                 else _number(self.unmodelled_wt_percent, 'unmodelled_wt_percent', 0, 100))
        if self.species_moles is not None:
            if ash or other:
                raise ValueError('For a species-defined feed, leave ash_wt_percent and '
                                 'unmodelled_wt_percent unset or zero.')
            return species_to_wt_percent(self.species_moles)
        if self.elements_wt_percent is None:
            raise ValueError('Give the feed as elements_wt_percent or species_moles.')
        analysis = {e: _number(v, f'feed {e} wt%', 0, 100) for e, v in self.elements_wt_percent.items()}
        if self.basis == 'dry':
            dry = analysis
        elif self.basis == 'dry_ash_free':
            if other:
                raise ValueError('A dry_ash_free analysis requires unmodelled_wt_percent = 0.')
            if not math.isclose(sum(analysis.values()), 100, abs_tol=0.05):
                raise ValueError('A dry_ash_free elemental analysis must sum to 100 %.')
            dry = {e: v * (1 - (ash or 0) / 100) for e, v in analysis.items()}
        else:
            raise ValueError("Feed basis must be 'dry' or 'dry_ash_free'.")
        if ash is not None or other is not None:
            total = sum(dry.values()) + (ash or 0) + (other or 0)
            if not math.isclose(total, 100, abs_tol=0.05):
                raise ValueError(f'Dry elements + ash + unmodelled = {total:.4g} %, not 100 %. '
                                 'No automatic renormalization is applied.')
        return dry


@dataclass
class Case:
    """One equilibrium calculation, or a sweep.

    temperature_C and pressure_bar may each be a number or a list; lists give a
    full grid (every temperature at every pressure). Alternatively ``points``
    lists explicit [temperature_C, pressure_bar] pairs (then set the other two
    to None).
    """
    temperature_C: float | list | None
    pressure_bar: float | list | None
    feed: Feed
    phases: dict
    name: str = 'case'
    description: str = ''
    # Atomic-weight overrides [g/mol]; change feed and balance arithmetic only.
    atomic_weights: dict = field(default_factory=dict)

    # Aqueous phase
    aqueous_model: str = 'electrolyte'  # 'electrolyte' (ions allowed) or 'molecular'
    aqueous_enabled: bool = True  # False removes the aqueous phase at every temperature
    aqueous_cutoff_C: float | None = None  # aqueous phase removed at T >= cutoff (a modelling rule)
    aqueous_water_floor: float = 1e-8  # minimum aqueous water [mol/mol dry-feed C]

    # Constraints. Carbon shares count carbon atoms, relative to ALL feed
    # carbon (including salt carbon); they are allocation rules, not kinetics.
    # carbon_bounds: {phase: {'exact': f} or {'min': f, 'max': f}} (fractions 0-1);
    # 'solid' means solid carbon only or all solids (solid_carbon_target).
    carbon_bounds: dict = field(default_factory=dict)
    solid_carbon_target: str = 'char'  # 'char' (solid carbon only) or 'all_solids'
    minimum_amounts: dict = field(default_factory=dict)  # {(phase, species): mol/mol C}
    # Old (format 1) carbon settings. They are moved into carbon_bounds when the
    # case is created (see legacy_carbon_bounds) and are None afterwards.
    gas_carbon_fraction: float | None = None  # exact fraction of feed C in gas
    gas_carbon_min_fraction: float | None = None  # at least this fraction in gas
    solid_carbon_fraction: float | None = None  # exact fraction in solids
    phase_carbon_fractions: dict | None = None  # {phase: exact fraction}

    # Starting point
    initialization: str = 'auto'
    # 'auto'   aqueous-water-rich start when an aqueous phase exists, else gas-water-rich
    # 'liquid' aqueous-water-rich start; plain feasible start without aqueous phase
    # 'vapor'  gas-water-rich start
    # 'lp'     any feasible start (zero LP objective)
    initial_gas_water_fraction: float = 1e-6  # cap on starting gas water for liquid starts
    initial_amounts: dict | None = None  # {(phase, species): mol/mol C}, replaces the LP start

    # Solver
    solver: str = 'ipopt'  # 'ipopt', 'slsqp', 'trust-constr', or 'slsqp-fd' (original)
    starts: int = 1  # starting points per point in total (initial, phase-rich/lean LP, random)
    distinct_composition_tol: float = 1e-4  # results differ above this |delta amount| [mol/mol C]
    distinct_gibbs_tol: float = 1.0  # ... or above this |delta G| [J/mol C]
    seed: int = 0
    stability_test: bool | None = None  # tangent-plane test of absent phases; default: starts > 1
    max_iterations: int = 5000
    energy_scale: float | None = None  # objective scale [J/mol C]; default max(1e4, |G_start|)
    fd_ftol: float = 1e-9  # 'slsqp-fd' only: tolerance on the scaled objective
    fd_eps: float = 1e-7  # 'slsqp-fd' only: finite-difference step
    solver_options: dict = field(default_factory=dict)  # backend-specific, e.g. ipopt_options
    # Old (format 1): extra random starts. Moved into ``starts`` when the case is
    # created (1 + 2 F + multistart, see starts_from_multistart); None afterwards.
    multistart: int | None = None

    # Sweeps
    # continuation: start each point from the previous solution. In a grid the
    # chain runs along temperature at each pressure; with explicit points it
    # follows the order of the list.
    continuation: bool = True
    refine: dict | None = None  # {start_C, stop_C, step_C}: extra temperatures in a grid
    points: list | None = None  # explicit [[T_C, P_bar], ...], replaces the grid

    # Data to plot next to the results: {file, quantity, label}. The CSV has a
    # T_C column and one column per 'phase:species'; quantity is a column of
    # the species table (e.g. 'mol_per_kg_dry_feed').
    comparison: dict | None = None

    # The feed as written in a format-2 case file (total mass, solids, ash,
    # elements as entered, or species amounts); ``feed`` is derived from it.
    feed_v2: dict | None = None

    def __post_init__(self):
        if self.feed is None and self.feed_v2 is not None:
            from .feed_v2 import feed_v2_to_model
            self.feed = feed_v2_to_model(self.feed_v2)
        elif isinstance(self.feed, dict):
            self.feed = Feed(**self.feed)
        self._move_legacy_settings()
        self.minimum_amounts = _pairs(self.minimum_amounts)
        if self.initial_amounts is not None:
            self.initial_amounts = _pairs(self.initial_amounts)
        self.phases = {phase: list(names) for phase, names in self.phases.items()}
        if self.points is not None:
            self.points = [list(p) for p in self.points]

    def _move_legacy_settings(self):
        """Old carbon settings -> carbon_bounds; multistart -> starts (the one place
        where format-1 settings are translated; the calculation reads only the new ones)."""
        legacy = {k: getattr(self, k) for k in LEGACY_CARBON_KEYS}
        if any(v not in (None, {}) for v in legacy.values()):
            bounds, self.solid_carbon_target = legacy_carbon_bounds(
                **legacy, solid_carbon_target=self.solid_carbon_target)
            self.carbon_bounds = merge_carbon_bounds(self.carbon_bounds, bounds)
        for key in LEGACY_CARBON_KEYS:
            setattr(self, key, None)
        self.carbon_bounds = normalize_carbon_bounds(self.carbon_bounds)
        if self.multistart is not None:
            if self.starts != 1:
                raise ValueError('Give either starts or the old multistart, not both.')
            try:
                temperatures = self.temperatures()
            except (ValueError, TypeError, KeyError):
                temperatures = None
            self.starts = starts_from_multistart(
                self.multistart, fluid_phase_count(self.phases, self.aqueous_enabled,
                                                   self.aqueous_cutoff_C, temperatures))
            self.multistart = None

    # ---------------------------------------------------------------- grid
    def temperatures(self):
        """Sorted temperatures [°C] of the case.

        In a grid this includes the refinement points and the aqueous cut-off
        temperature when it lies inside the temperature range.
        """
        if self.points is not None:
            return sorted(set(t for t, _ in self._explicit_points()))
        if self.temperature_C is None:
            raise ValueError('Give temperature_C (or points).')
        values = self.temperature_C if isinstance(self.temperature_C, (list, tuple)) else [self.temperature_C]
        grid = [_number(t, 'temperature_C', -273.149999) for t in values]
        if not grid:
            raise ValueError('Give at least one temperature.')
        grid = refined_temperatures(grid, self.refine)
        if self.aqueous_cutoff_C is not None and min(grid) < float(self.aqueous_cutoff_C) < max(grid):
            grid.append(float(self.aqueous_cutoff_C))
        return sorted(set(round(float(t), 8) for t in grid))

    def pressures(self):
        """Sorted pressures [bar] of the case."""
        if self.points is not None:
            return sorted(set(p for _, p in self._explicit_points()))
        if self.pressure_bar is None:
            raise ValueError('Give pressure_bar (or points).')
        values = self.pressure_bar if isinstance(self.pressure_bar, (list, tuple)) else [self.pressure_bar]
        if not values:
            raise ValueError('Give at least one pressure.')
        return sorted(set(round(_number(p, 'pressure_bar', 1e-300), 8) for p in values))

    def _explicit_points(self):
        out = []
        for item in self.points:
            if not isinstance(item, (list, tuple)) or len(item) != 2:
                raise ValueError(f'Each entry of points must be [temperature_C, pressure_bar], not {item!r}.')
            point = (round(_number(item[0], 'points temperature_C', -273.149999), 8),
                     round(_number(item[1], 'points pressure_bar', 1e-300), 8))
            if point not in out:
                out.append(point)
        if not out:
            raise ValueError('points is empty.')
        return out

    def grid(self):
        """All (temperature_C, pressure_bar) points: explicit points in their order,
        else every temperature at every pressure (by pressure, then temperature)."""
        if self.points is not None:
            return self._explicit_points()
        temperatures = self.temperatures()
        return [(T, P) for P in self.pressures() for T in temperatures]

    def chains(self):
        """Groups of points that must be solved in order (continuation).

        Without continuation every point is its own group. With continuation a
        grid gives one chain per pressure (rising temperature) and explicit
        points one chain in list order.
        """
        grid = self.grid()
        if not self.continuation:
            return [[point] for point in grid]
        if self.points is not None:
            return [grid]
        chains = {}
        for T, P in grid:
            chains.setdefault(P, []).append((T, P))
        return list(chains.values())

    def is_single_point(self):
        return len(self.grid()) == 1

    def at(self, temperature_C, pressure_bar=None):
        """Copy of this case at one temperature (and pressure; optional when the
        case has a single pressure)."""
        if pressure_bar is None:
            pressures = self.pressures()
            if len(pressures) != 1:
                raise ValueError('This case has several pressures; give pressure_bar too.')
            pressure_bar = pressures[0]
        point = copy.deepcopy(self)
        point.temperature_C = float(temperature_C)
        point.pressure_bar = float(pressure_bar)
        point.points = None
        point.refine = None
        return point

    def replace(self, **changes):
        point = copy.deepcopy(self)
        for key, value in changes.items():
            if not hasattr(point, key):
                raise AttributeError(f'Case has no field {key!r}')
            setattr(point, key, value)
        if 'feed' in changes and 'feed_v2' not in changes:
            point.feed_v2 = None  # the written feed no longer describes the new one
        elif 'feed_v2' in changes and 'feed' not in changes and point.feed_v2 is not None:
            from .feed_v2 import feed_v2_to_model
            point.feed = feed_v2_to_model(point.feed_v2, salts=self.feed.salts if self.feed else None)
        point.__post_init__()
        return point

    def validate(self):
        """Checks that need no databases. Raises ValueError with a readable message."""
        unknown = set(self.phases) - set(PHASES)
        if unknown:
            raise ValueError(f'Unknown phase(s) {sorted(unknown)}; use {list(PHASES)}.')
        if self.solver not in SOLVERS:
            raise ValueError(f'solver must be one of {SOLVERS}, not {self.solver!r}.')
        if self.initialization not in INITIALIZATIONS:
            raise ValueError(f'initialization must be one of {INITIALIZATIONS}.')
        if self.aqueous_model not in ('electrolyte', 'molecular'):
            raise ValueError("aqueous_model must be 'electrolyte' or 'molecular'.")
        if self.solid_carbon_target not in ('char', 'all_solids'):
            raise ValueError("solid_carbon_target must be 'char' or 'all_solids'.")
        if self.points is not None and self.refine:
            raise ValueError('refine applies to a temperature grid; it cannot be used with points.')
        check_carbon_bounds(self.carbon_bounds)
        if isinstance(self.starts, bool) or not isinstance(self.starts, (int, float)) \
                or self.starts != int(self.starts) or self.starts < 1:
            raise ValueError(f'starts (initial guesses) must be a whole number of at least 1, '
                             f'not {self.starts!r}.')
        self.starts = int(self.starts)
        for name in ('distinct_composition_tol', 'distinct_gibbs_tol'):
            _number(getattr(self, name), name, 0)
        self.grid()
        return self

    # ---------------------------------------------------------------- files
    def to_toml(self):
        return dump_case(self)

    def save(self, path):
        Path(path).write_text(self.to_toml(), encoding='utf-8')


def _number(value, name, low=None, high=None):
    try:
        value = float(value)
    except (TypeError, ValueError):
        raise ValueError(f'{name} must be a number, not {value!r}.') from None
    if not math.isfinite(value) or (low is not None and value < low) or (high is not None and value > high):
        raise ValueError(f'Invalid {name}: {value}')
    return value


def _pairs(mapping):
    """{(phase, species): v} from either that form or nested {phase: {species: v}}."""
    out = {}
    for key, value in (mapping or {}).items():
        if isinstance(key, tuple):
            out[key] = float(value)
        else:
            for species, amount in value.items():
                out[(key, species)] = float(amount)
    return out


def _nested(pairs):
    out = {}
    for (phase, species), value in pairs.items():
        out.setdefault(phase, {})[species] = value
    return out


# ---------------------------------------------------------------- carbon bounds
PHASE_LABELS = {'gas': 'Gas', 'aqueous': 'Aqueous', 'organic_liquid': 'Organic liquid', 'solid': 'Solid'}
FLUID_PHASES = ('gas', 'aqueous', 'organic_liquid')
BOUND_KINDS = ('exact', 'min', 'max')
LEGACY_CARBON_KEYS = ('gas_carbon_fraction', 'gas_carbon_min_fraction', 'solid_carbon_fraction',
                      'phase_carbon_fractions')


def _label(phase):
    return PHASE_LABELS.get(phase, str(phase))


def legacy_carbon_bounds(gas_carbon_fraction=None, gas_carbon_min_fraction=None,
                         solid_carbon_fraction=None, phase_carbon_fractions=None,
                         solid_carbon_target='char'):
    """(carbon_bounds, solid_carbon_target) for the old carbon settings.

    gas_carbon_fraction -> gas.exact, gas_carbon_min_fraction -> gas.min,
    solid_carbon_fraction -> solid.exact, phase_carbon_fractions[p] -> p.exact.
    The old phase_carbon_fractions['solid'] counted all solids, so it sets
    solid_carbon_target = 'all_solids' (it cannot be combined with
    solid_carbon_fraction, which counted char only by default).
    """
    bounds = {}
    if gas_carbon_fraction is not None and gas_carbon_min_fraction is not None:
        raise ValueError('Use either gas_carbon_fraction (equality) or '
                         'gas_carbon_min_fraction (lower bound), not both')
    if gas_carbon_fraction is not None:
        bounds['gas'] = {'exact': gas_carbon_fraction}
    if gas_carbon_min_fraction is not None:
        bounds['gas'] = {'min': gas_carbon_min_fraction}
    if solid_carbon_fraction is not None:
        bounds['solid'] = {'exact': solid_carbon_fraction}
    for phase, fraction in (phase_carbon_fractions or {}).items():
        if phase == 'gas' and gas_carbon_fraction is not None:
            if not math.isclose(float(fraction), float(gas_carbon_fraction), rel_tol=0, abs_tol=1e-10):
                raise ValueError('Conflicting gas-carbon constraints')
            continue
        if phase == 'gas' and gas_carbon_min_fraction is not None:
            raise ValueError('Conflicting gas-carbon constraints')
        if phase == 'solid':
            if solid_carbon_fraction is not None:
                raise ValueError('phase_carbon_fractions.solid (all solids) cannot be combined with '
                                 'solid_carbon_fraction; use [constraints.carbon] solid and '
                                 'solid_carbon_target.')
            solid_carbon_target = 'all_solids'
        bounds[phase] = {'exact': fraction}
    return bounds, solid_carbon_target


def merge_carbon_bounds(bounds, extra):
    out = {phase: dict(b) if isinstance(b, dict) else b for phase, b in (bounds or {}).items()}
    for phase, bound in extra.items():
        if phase in out and out[phase] != bound:
            raise ValueError(f'{_label(phase)}: the carbon share is given twice (old and new '
                             'keys) with different values.')
        out[phase] = dict(bound)
    return out


def normalize_carbon_bounds(bounds):
    """Ordered copy: gas, solid, then the other phases as given; numbers as floats."""
    bounds = dict(bounds or {})
    order = [p for p in ('gas', 'solid') if p in bounds] + [p for p in bounds if p not in ('gas', 'solid')]
    out = {}
    for phase in order:
        bound = bounds[phase]
        if not isinstance(bound, dict):
            out[phase] = bound
            continue
        keys = [k for k in BOUND_KINDS if k in bound] + [k for k in bound if k not in BOUND_KINDS]
        item = {}
        for key in keys:
            value = bound[key]
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                value = float(value)
            item[key] = value
        out[phase] = item
    return out


def check_carbon_bounds(bounds):
    """Raise ValueError for the first malformed bound (needs no species data)."""
    for phase, bound in (bounds or {}).items():
        if phase not in PHASES:
            raise ValueError(f'Unknown phase {phase!r} in [constraints.carbon]; use '
                             f'{", ".join(PHASES)}.')
        label = _label(phase)
        if not isinstance(bound, dict) or not bound or set(bound) - set(BOUND_KINDS):
            raise ValueError(f'{label} carbon: use exact, or min and/or max '
                             f'(e.g. {phase} = {{ min = 0.9 }}).')
        if 'exact' in bound and ('min' in bound or 'max' in bound):
            raise ValueError(f'{label} carbon: exact cannot be combined with min or max.')
        for kind, value in bound.items():
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) \
                    or not 0 <= value <= 1:
                raise ValueError(f'{label} carbon: {kind} must be a fraction between 0 and 1 '
                                 f'(0 % to 100 %), not {value!r}.')
        if 'min' in bound and 'max' in bound and bound['min'] > bound['max']:
            raise ValueError(f'{label}: the minimum is larger than the maximum.')


def carbon_bound_messages(case, carbon_solids=None):
    """Readable problems with the carbon bounds of ``case`` (empty list: none found).

    carbon_solids  optional {solid species name: True if it contains carbon}; used
                   to decide whether "char only" covers all solid carbon.
    """
    messages = []
    bounds = case.carbon_bounds
    try:
        check_carbon_bounds(bounds)
    except ValueError as exc:
        return [str(exc)]
    present = {p for p in PHASES if case.phases.get(p)}
    if case.aqueous_enabled is False:
        present.discard('aqueous')
    for phase in bounds:
        if phase not in present:
            messages.append(f'{_label(phase)} has a carbon bound but no species.')
    fixed = sum(b.get('exact', 0.0) + b.get('min', 0.0) for b in bounds.values())
    if fixed > 1 + 1e-10:
        messages.append('The carbon shares add up to more than 100 %.')
    elif present and all('exact' in bounds.get(p, {}) for p in present):
        solid_complete = 'solid' not in present or case.solid_carbon_target == 'all_solids' or (
            carbon_solids is not None
            and not any(carbon_solids.get(s, True) for s in case.phases['solid'] if s != 'carbon'))
        total = sum(bounds[p]['exact'] for p in present)
        if solid_complete and not math.isclose(total, 1.0, rel_tol=0, abs_tol=1e-10):
            messages.append('With every phase fixed, the carbon shares must add up to 100 %.')
    return messages


def effective_carbon_bounds(bounds, aqueous_disabled):
    """(bounds, dropped) at a point: without an aqueous phase the aqueous bound is
    dropped; dropped = ['aqueous carbon'] when that bound asked for carbon there."""
    if not aqueous_disabled or 'aqueous' not in bounds:
        return dict(bounds), []
    bound = bounds['aqueous']
    out = {p: b for p, b in bounds.items() if p != 'aqueous'}
    asked = bound.get('exact', 0.0) > 0 or bound.get('min', 0.0) > 0
    return out, ['aqueous carbon'] if asked else []


# ---------------------------------------------------------------- starting points
def fluid_phase_count(phases, aqueous_enabled=True, aqueous_cutoff_C=None, temperatures=None):
    """Fluid phases (gas, aqueous, organic liquid) with at least one species.

    The aqueous phase does not count when it is disabled, or removed by the
    cut-off at every temperature of the case.
    """
    count = 0
    for phase in FLUID_PHASES:
        if not (phases or {}).get(phase):
            continue
        if phase == 'aqueous':
            if aqueous_enabled is False:
                continue
            if aqueous_cutoff_C is not None and temperatures and \
                    min(temperatures) >= float(aqueous_cutoff_C):
                continue
        count += 1
    return count


def starts_from_multistart(multistart, fluid_phases):
    """Total starting points of the old multistart = k: 1 + 2 F + k (k = 0: 1).

    The old code made the initial start, a phase-rich and a phase-lean LP start
    for each of the F fluid phases, and k random starts.
    """
    k = int(multistart or 0)
    if k < 0:
        raise ValueError('multistart cannot be negative.')
    return 1 + 2 * int(fluid_phases) + k if k > 0 else 1


def refined_temperatures(grid, refine):
    """``grid`` plus the extra temperatures of the old ``refine`` setting
    ({start_C, stop_C, step_C}, limited to the range of the grid)."""
    grid = list(grid)
    if refine:
        step = _number(refine['step_C'], 'refine step_C', 1e-12)
        a = max(min(grid), float(refine['start_C']))
        b = min(max(grid), float(refine['stop_C']))
        if b >= a:
            grid += [float(t) for t in np.arange(a, b + step * 1e-8, step)]
    return grid


# ---------------------------------------------------------------- old (format 1) files
class LegacyConversionError(ValueError):
    """An old case file cannot be written in format 2; the message says why."""


NEW_KEYS = {'feed': ('total_kg', 'solids', 'ash', 'species'), 'constraints': ('carbon',),
            'solver': ('starts', 'distinct_composition_tol', 'distinct_gibbs_tol')}


def is_v2(data):
    """True for case-file data in format 2 (format = 2, or any key only format 2 has)."""
    if data.get('format') == 2:
        return True
    return any(isinstance(data.get(section), dict) and key in data[section]
               for section, keys in NEW_KEYS.items() for key in keys)


def legacy_to_v2(data, molar_mass=None):
    """Case-file data in format 2 for old (format-1) data, giving the same calculation.

    Raises LegacyConversionError when the file cannot be represented (its message
    is meant for the user): unmodelled_wt_percent > 0, or a species-moles feed
    without ``molar_mass(name) -> g/mol``. Data already in format 2 is returned
    as a copy. See htg.feed_v2.legacy_feed_to_v2 for the feed rules.
    """
    from .feed_v2 import legacy_feed_to_v2
    data = copy.deepcopy(data)
    if is_v2(data):
        return data
    out = {'format': 2}
    out.update(data)
    if 'feed' in out:
        out['feed'] = legacy_feed_to_v2(out['feed'], molar_mass)
    constraints = out.get('constraints')
    if isinstance(constraints, dict):
        legacy = {key: constraints.pop(key) for key in LEGACY_CARBON_KEYS if key in constraints}
        if legacy:
            bounds, target = legacy_carbon_bounds(
                **legacy, solid_carbon_target=constraints.get('solid_carbon_target', 'char'))
            if target != constraints.get('solid_carbon_target', 'char'):
                constraints['solid_carbon_target'] = target
            if bounds:
                constraints['carbon'] = normalize_carbon_bounds(bounds)
    solver = out.get('solver')
    if isinstance(solver, dict) and 'multistart' in solver:
        conditions = out.get('conditions', {})
        aqueous = out.get('aqueous', {})
        temperatures = None
        if 'temperature_C' in conditions:
            values = conditions['temperature_C']
            temperatures = refined_temperatures(values if isinstance(values, list) else [values],
                                                (out.get('sweep') or {}).get('refine'))
        elif 'points' in conditions:
            temperatures = [t for t, _ in conditions['points']]
        count = fluid_phase_count(out.get('phases', {}), aqueous.get('enabled', True),
                                  aqueous.get('cutoff_C'), temperatures)
        solver['starts'] = starts_from_multistart(solver.pop('multistart'), count)
    sweep = out.get('sweep')
    if isinstance(sweep, dict) and 'refine' in sweep:
        refine = sweep.pop('refine')
        conditions = out.get('conditions', {})
        if refine and 'temperature_C' in conditions and 'points' not in conditions:
            values = conditions['temperature_C']
            values = values if isinstance(values, list) else [values]
            merged = sorted(set(round(float(t), 8) for t in refined_temperatures(values, refine)))
            conditions['temperature_C'] = [int(t) if t.is_integer() else t for t in merged]
        if not sweep:
            out.pop('sweep')
    return out


# ---------------------------------------------------------------- TOML case files
# section -> {key in file: Case field}. Old (format 1) keys are still read.
SECTIONS = {
    'conditions': {'temperature_C': 'temperature_C', 'pressure_bar': 'pressure_bar', 'points': 'points'},
    'aqueous': {'model': 'aqueous_model', 'enabled': 'aqueous_enabled', 'cutoff_C': 'aqueous_cutoff_C',
                'water_floor': 'aqueous_water_floor'},
    'constraints': {'carbon': 'carbon_bounds', 'solid_carbon_target': 'solid_carbon_target',
                    'minimum_amounts': 'minimum_amounts',
                    **{k: k for k in LEGACY_CARBON_KEYS}},
    'initial_guess': {'initialization': 'initialization',
                      'initial_gas_water_fraction': 'initial_gas_water_fraction',
                      'amounts': 'initial_amounts'},
    'solver': {'method': 'solver', 'starts': 'starts',
               'distinct_composition_tol': 'distinct_composition_tol',
               'distinct_gibbs_tol': 'distinct_gibbs_tol', 'seed': 'seed',
               'stability_test': 'stability_test', 'max_iterations': 'max_iterations',
               'energy_scale': 'energy_scale', 'fd_ftol': 'fd_ftol', 'fd_eps': 'fd_eps',
               'options': 'solver_options', 'multistart': 'multistart'},
    'sweep': {'continuation': 'continuation', 'refine': 'refine'},
    'comparison': {'file': 'file', 'quantity': 'quantity', 'label': 'label'},
}
TOP_LEVEL = ('name', 'description', 'atomic_weights')
FORMAT = 2
LEGACY_FEED_ONLY = ('basis', 'ash_wt_percent', 'unmodelled_wt_percent', 'species_moles',
                    'moisture_wt_percent', 'water_g_per_g_dry')


def _feed_from_dict(feed, version):
    """(Feed, feed_v2) for a [feed] table of either format."""
    from .feed_v2 import feed_v2_to_model, is_v2_feed
    if not isinstance(feed, dict):
        raise ValueError('[feed] must be a table.')
    if is_v2_feed(feed):
        mixed = [k for k in LEGACY_FEED_ONLY if k in feed]
        if mixed:
            raise ValueError(f'[feed] mixes the new keys (total_kg, solids, ash, species) with old '
                             f'ones {mixed}; remove the old ones.')
        written = {k: copy.deepcopy(v) for k, v in feed.items() if k != 'salts'}
        return feed_v2_to_model(feed), written
    if version == FORMAT and not any(k in feed for k in LEGACY_FEED_ONLY):
        raise ValueError('[feed] needs total_kg and solids (and elements_wt_percent), '
                         'or a [feed.species] table.')
    known = {f.name for f in fields(Feed)} - {'species_amounts'}
    unknown = set(feed) - known
    if unknown:
        raise ValueError(f'Unknown [feed] key(s) {sorted(unknown)}; allowed: {sorted(known)}.')
    return Feed(**feed), None


def case_from_dict(data, name=None):
    """Case from case-file data, format 2 or the old format (or a mixture)."""
    data = copy.deepcopy(data)
    version = data.pop('format', None)
    if version not in (None, 1, FORMAT):
        raise ValueError(f'format = {version!r} is not known; this version reads format 2 '
                         '(and old files without format).')
    kwargs = {}
    for key in TOP_LEVEL:
        if key in data:
            kwargs[key] = data.pop(key)
    if 'feed' not in data:
        raise ValueError('The case needs a [feed] section.')
    kwargs['feed'], kwargs['feed_v2'] = _feed_from_dict(data.pop('feed'), version)
    kwargs['phases'] = data.pop('phases', {})
    for section, keys in SECTIONS.items():
        values = data.pop(section, {})
        if section == 'comparison':
            if values:
                unknown = set(values) - set(keys)
                if unknown or 'file' not in values:
                    raise ValueError('[comparison] needs file (and optionally quantity, label).')
                kwargs['comparison'] = values
            continue
        unknown = set(values) - set(keys)
        if unknown:
            raise ValueError(f'Unknown [{section}] key(s) {sorted(unknown)}; allowed: {sorted(keys)}.')
        for key, value in values.items():
            kwargs[keys[key]] = value
    if data:
        raise ValueError(f'Unknown section(s) or key(s): {sorted(data)}.')
    if 'points' in kwargs:
        if 'temperature_C' in kwargs or 'pressure_bar' in kwargs:
            raise ValueError('[conditions] has points; leave out temperature_C and pressure_bar.')
        kwargs['temperature_C'] = kwargs['pressure_bar'] = None
    else:
        for required in ('temperature_C', 'pressure_bar'):
            if required not in kwargs:
                raise ValueError(f'[conditions] needs {required} (or points).')
    if name and 'name' not in kwargs:
        kwargs['name'] = name
    return Case(**kwargs).validate()


def load_case(path):
    """Read a TOML case file."""
    path = Path(path)
    with open(path, 'rb') as f:
        case = case_from_dict(tomllib.load(f), name=path.stem)
    if case.comparison:
        case.comparison['file'] = str((path.parent / case.comparison['file']).resolve())
    return case


def _feed_for_file(case):
    """The [feed] table to write: format 2 when possible, else the old keys."""
    from .feed_v2 import legacy_feed_to_v2
    if case.feed_v2 is not None:
        return copy.deepcopy(case.feed_v2)
    feed = asdict(case.feed)
    feed.pop('salts')
    if feed.get('species_amounts') is not None:
        return {'species': feed['species_amounts']}
    feed.pop('species_amounts')
    try:
        written = legacy_feed_to_v2(feed)
    except LegacyConversionError:
        # species moles (needs molar masses) or unmodelled mass: keep the old keys
        defaults = {f.name: f.default for f in fields(Feed) if f.name not in ('salts', 'species_amounts')}
        return {k: v for k, v in feed.items() if v != defaults[k]}
    written.pop('salts', None)
    return written


def dump_case(case):
    """TOML text (format 2) for ``case``: values that differ from the defaults, plus required ones."""
    defaults = {f.name: f.default for f in fields(Case) if f.default is not MISSING}
    lines = [f'format = {FORMAT}', f'name = {_toml(case.name)}']
    if case.description:
        lines.append(f'description = {_toml(case.description)}')
    if case.atomic_weights:
        lines.append(f'atomic_weights = {_toml(case.atomic_weights)}')

    def section(title, items):
        items = [(k, v) for k, v in items if v is not None]
        if items:
            lines.append(f'\n[{title}]')
            lines.extend(f'{_key(k)} = {_toml(v)}' for k, v in items)

    temperatures = case.temperature_C
    if case.refine and case.points is None and temperatures is not None:
        values = temperatures if isinstance(temperatures, (list, tuple)) else [temperatures]
        merged = sorted(set(round(float(t), 8) for t in refined_temperatures(values, case.refine)))
        temperatures = [int(t) if t.is_integer() else t for t in merged]
    section('conditions', [('temperature_C', temperatures), ('pressure_bar', case.pressure_bar),
                           ('points', case.points)])
    feed = _feed_for_file(case)
    species = feed.pop('species', None)
    section('feed', list(feed.items()))
    if species is not None:
        section('feed.species', list(species.items()))
    for salt_name, salt in case.feed.salts.items():
        lines.append(f'\n[feed.salts.{_key(salt_name)}]')
        lines.append(f'mol_per_g_dry_feed = {_toml(salt.mol_per_g_dry_feed)}')
        lines.append(f'elements = {_toml(salt.elements)}')
    section('phases', [(p, case.phases[p]) for p in PHASES if p in case.phases])
    for title, keys in SECTIONS.items():
        if title == 'conditions':
            continue
        if title == 'comparison':
            if case.comparison:
                section(title, list(case.comparison.items()))
            continue
        items = []
        for key, attr in keys.items():
            if attr in ('refine', 'carbon_bounds'):
                continue
            value = getattr(case, attr)
            if attr in ('minimum_amounts', 'initial_amounts') and value:
                value = _nested(value)
            if attr in defaults and value == defaults[attr]:
                continue
            if value in ({}, None):
                continue
            items.append((key, value))
        section(title, items)
        if title == 'constraints' and case.carbon_bounds:
            section('constraints.carbon', list(case.carbon_bounds.items()))
    return '\n'.join(lines) + '\n'


def _key(text):
    return text if text.replace('_', '').replace('-', '').isalnum() else _toml(text)


def _toml(value):
    if isinstance(value, bool):
        return 'true' if value else 'false'
    if isinstance(value, (int, float)):
        return repr(value)
    if isinstance(value, str):
        return '"' + value.replace('\\', '\\\\').replace('"', '\\"') + '"'
    if isinstance(value, (list, tuple)):
        return '[' + ', '.join(_toml(v) for v in value) + ']'
    if isinstance(value, dict):
        return '{ ' + ', '.join(f'{_key(str(k))} = {_toml(v)}' for k, v in value.items()) + ' }'
    raise TypeError(f'Cannot write {value!r} to TOML')
