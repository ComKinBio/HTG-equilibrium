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
from dataclasses import asdict, dataclass, field, fields
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
    """
    elements_wt_percent: dict | None = None
    basis: str = 'dry'
    ash_wt_percent: float | None = None
    unmodelled_wt_percent: float | None = None
    species_moles: dict | None = None
    moisture_wt_percent: float = 0.0
    water_g_per_g_dry: float = 0.0
    salts: dict = field(default_factory=dict)

    def __post_init__(self):
        self.salts = {name: s if isinstance(s, Salt) else Salt(**s) for name, s in self.salts.items()}

    def total_water_g_per_g_dry(self):
        moisture = _number(self.moisture_wt_percent, 'feed moisture_wt_percent', 0)
        if moisture >= 100:
            raise ValueError('Feed moisture must be below 100 wt%.')
        return _number(self.water_g_per_g_dry, 'feed water_g_per_g_dry', 0) + moisture / (100 - moisture)

    def dry_element_wt_percent(self, species_to_wt_percent):
        """Dry-basis elemental wt% used by the model.

        ``species_to_wt_percent`` converts ``species_moles`` (it needs the registry).
        """
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

    # Constraints. Carbon fractions count carbon atoms, relative to ALL feed
    # carbon (including salt carbon); they are allocation rules, not kinetics.
    gas_carbon_fraction: float | None = None  # exact fraction of feed C in gas
    gas_carbon_min_fraction: float | None = None  # at least this fraction in gas
    solid_carbon_fraction: float | None = None  # exact fraction in solids
    solid_carbon_target: str = 'char'  # 'char' (solid carbon only) or 'all_solids'
    phase_carbon_fractions: dict = field(default_factory=dict)  # {phase: exact fraction}
    minimum_amounts: dict = field(default_factory=dict)  # {(phase, species): mol/mol C}

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
    multistart: int = 0  # extra random feasible starts (plus phase-rich starts) when > 0
    seed: int = 0
    stability_test: bool | None = None  # tangent-plane test of absent phases; default: with multistart
    max_iterations: int = 5000
    energy_scale: float | None = None  # objective scale [J/mol C]; default max(1e4, |G_start|)
    fd_ftol: float = 1e-9  # 'slsqp-fd' only: tolerance on the scaled objective
    fd_eps: float = 1e-7  # 'slsqp-fd' only: finite-difference step
    solver_options: dict = field(default_factory=dict)  # backend-specific, e.g. ipopt_options

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

    def __post_init__(self):
        if isinstance(self.feed, dict):
            self.feed = Feed(**self.feed)
        self.minimum_amounts = _pairs(self.minimum_amounts)
        if self.initial_amounts is not None:
            self.initial_amounts = _pairs(self.initial_amounts)
        self.phases = {phase: list(names) for phase, names in self.phases.items()}
        if self.points is not None:
            self.points = [list(p) for p in self.points]

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
        if self.refine:
            step = _number(self.refine['step_C'], 'refine step_C', 1e-12)
            a = max(min(grid), float(self.refine['start_C']))
            b = min(max(grid), float(self.refine['stop_C']))
            if b >= a:
                grid += [float(t) for t in np.arange(a, b + step * 1e-8, step)]
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


# ---------------------------------------------------------------- TOML case files
# section -> {key in file: Case field}
SECTIONS = {
    'conditions': {'temperature_C': 'temperature_C', 'pressure_bar': 'pressure_bar', 'points': 'points'},
    'aqueous': {'model': 'aqueous_model', 'enabled': 'aqueous_enabled', 'cutoff_C': 'aqueous_cutoff_C',
                'water_floor': 'aqueous_water_floor'},
    'constraints': {k: k for k in ('gas_carbon_fraction', 'gas_carbon_min_fraction',
                                   'solid_carbon_fraction', 'solid_carbon_target',
                                   'phase_carbon_fractions', 'minimum_amounts')},
    'initial_guess': {'initialization': 'initialization',
                      'initial_gas_water_fraction': 'initial_gas_water_fraction',
                      'amounts': 'initial_amounts'},
    'solver': {'method': 'solver', 'multistart': 'multistart', 'seed': 'seed',
               'stability_test': 'stability_test', 'max_iterations': 'max_iterations',
               'energy_scale': 'energy_scale', 'fd_ftol': 'fd_ftol', 'fd_eps': 'fd_eps',
               'options': 'solver_options'},
    'sweep': {'continuation': 'continuation', 'refine': 'refine'},
    'comparison': {'file': 'file', 'quantity': 'quantity', 'label': 'label'},
}
TOP_LEVEL = ('name', 'description', 'atomic_weights')


def case_from_dict(data, name=None):
    data = copy.deepcopy(data)
    kwargs = {}
    for key in TOP_LEVEL:
        if key in data:
            kwargs[key] = data.pop(key)
    if 'feed' not in data:
        raise ValueError('The case needs a [feed] section.')
    feed = data.pop('feed')
    known = {f.name for f in fields(Feed)}
    unknown = set(feed) - known
    if unknown:
        raise ValueError(f'Unknown [feed] key(s) {sorted(unknown)}; allowed: {sorted(known)}.')
    kwargs['feed'] = Feed(**feed)
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


def dump_case(case):
    """TOML text for ``case`` (only values that differ from the defaults, plus required ones)."""
    defaults = {f.name: f.default for f in fields(Case) if f.default is not f.default_factory}
    data = asdict(case)
    lines = [f'name = {_toml(case.name)}']
    if case.description:
        lines.append(f'description = {_toml(case.description)}')
    if case.atomic_weights:
        lines.append(f'atomic_weights = {_toml(case.atomic_weights)}')

    def section(title, items):
        items = [(k, v) for k, v in items if v is not None]
        if items:
            lines.append(f'\n[{title}]')
            lines.extend(f'{k} = {_toml(v)}' for k, v in items)

    section('conditions', [('temperature_C', case.temperature_C), ('pressure_bar', case.pressure_bar),
                           ('points', case.points)])
    feed = {k: v for k, v in data['feed'].items() if k != 'salts'}
    feed_defaults = {f.name: f.default for f in fields(Feed) if f.name != 'salts'}
    section('feed', [(k, v) for k, v in feed.items() if v != feed_defaults[k]])
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
            value = getattr(case, attr)
            if attr in ('minimum_amounts', 'initial_amounts') and value:
                value = _nested(value)
            if attr in defaults and value == defaults[attr]:
                continue
            if value in ({}, None):
                continue
            items.append((key, value))
        section(title, items)
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
