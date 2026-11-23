"""Build PlotSpecs for the Results tab from result tables.

Every builder takes ``tables`` (any object with ``summary``, ``species`` and
``phases`` DataFrames, e.g. htg.gui.runs.RunTables) and PlotOptions, and returns a
htg.plotspec.PlotSpec. Values of temperature and pressure in the options are in
case-file units (°C, bar); axis values are converted for display with the
converter callables in the options (identity by default), so this module does not
depend on any unit module.
"""
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import quantities as q
from .plotspec import NEUTRAL, Heatmap, Panel, PlotSpec, Series, VLine, assign_colours

PLOTS = {'total': 'Amount, all phases together', 'phases': 'Amount in each phase',
         'gas': 'Gas composition', 'cge': 'Carbon gasification efficiency',
         'carbon': 'Carbon in each phase', 'map': 'T–P map'}
AMOUNT_PLOTS = ('total', 'phases')
SPECIES_PLOTS = ('total', 'phases', 'gas')
PANEL_ORDER = ('gas', 'aqueous', 'solid', 'organic_liquid')
PHASE_TITLES = {'gas': 'Gas', 'aqueous': 'Aqueous', 'solid': 'Solid', 'organic_liquid': 'Organic liquid'}
PHASE_COLOURS = {'gas': 0, 'aqueous': 1, 'organic_liquid': 2, 'solid': 3}
LARGEST = 8
LOG_FLOOR = 1e-8  # values at or below this are not drawn on log axes
SHOWN = 1e-6  # a phase line in "Carbon in each phase" needs more than this [%] (noise ~1e-7)
NO_POINTS = 'No points solved yet.'


def number(value):
    """Up to 6 significant figures, no trailing zeros."""
    return f'{float(value):.6g}'


def _identity(value):
    return value


@dataclass
class Measured:
    """Measured (reference) data: wide table T_C [, P_bar] + 'phase:species' columns."""
    table: pd.DataFrame
    quantity: str = 'mol_per_kg_dry_feed'
    label: str = 'Measured'


@dataclass
class PlotOptions:
    """Choices of the Results tab settings panel.

    x_axis        'T' (temperature) or 'P' (pressure)
    fixed         value of the other quantity in case units (°C or bar); None = first value
    amount        'mol' or 'g'; basis 'total', 'dry' or 'daf'; s and a as in htg.quantities
    species       species to draw (None = the largest 8); water is drawn only when
                  include_water is set
    composition   'mol' or 'mass' (gas composition); exclude_water for dry gas
    map_item      item of the T–P map (see map_items)
    cutoff_C      aqueous cut-off [°C] of the case; marked on temperature axes
    measured      optional Measured data for the amount plots
    colours       {species: palette index}; None = run_colours(species table)
    t_unit, p_unit, t_convert, p_convert, number_format
                  display units: labels, converters from °C / bar, and the number format
    """
    x_axis: str = 'T'
    fixed: float | None = None
    amount: str = 'mol'
    basis: str = 'dry'
    s: float | None = None
    a: float = 0.0
    species: list | None = None
    include_water: bool = False
    log: bool = False
    composition: str = 'mol'
    exclude_water: bool = False
    map_item: str = 'gas_carbon_percent'
    cutoff_C: float | None = None
    measured: Measured | None = None
    colours: dict | None = None
    t_unit: str = '°C'
    p_unit: str = 'bar'
    t_convert: Callable = field(default=_identity)
    p_convert: Callable = field(default=_identity)
    number_format: Callable = field(default=number)

    @property
    def axis(self):
        return 'T_C' if self.x_axis == 'T' else 'P_bar'

    @property
    def other(self):
        return 'P_bar' if self.x_axis == 'T' else 'T_C'

    def convert(self, column, value):
        return (self.t_convert if column == 'T_C' else self.p_convert)(float(value))

    def unit(self, column):
        return self.t_unit if column == 'T_C' else self.p_unit

    def text(self, column, value):
        """'240 bar' in the display unit."""
        return f'{self.number_format(self.convert(column, value))} {self.unit(column)}'

    def axis_label(self, column):
        return f"{'Temperature' if column == 'T_C' else 'Pressure'} [{self.unit(column)}]"


# --- helpers for the settings panel -------------------------------------------------------------

def run_colours(species):
    """{species: palette index} for one run, the same in every plot and panel.

    Species are ordered by their largest total amount (all phases, mol / mol C) so the
    main species get the first, unrepeated colours; water comes last.
    """
    if species is None or species.empty:
        return {}
    total = species.groupby(['T_C', 'P_bar', 'species'], sort=False).mol_per_mol_C.sum()
    peak = total.groupby(level='species', sort=False).max().sort_values(ascending=False, kind='stable')
    names = [s for s in peak.index if s != q.WATER]
    if q.WATER in peak.index:
        names.append(q.WATER)
    return assign_colours(names)


def species_choices(species):
    """Species with a nonzero amount somewhere in the run, in colour order."""
    colours = run_colours(species)
    present = set(species.loc[species.mol_per_mol_C.astype(float) > q.TINY, 'species'])
    return [s for s in colours if s in present]


def map_available(summary):
    """The T–P map needs at least 2 temperatures and 2 pressures."""
    return summary['T_C'].nunique() >= 2 and summary['P_bar'].nunique() >= 2 if len(summary) else False


def fixed_choices(tables, options):
    """Values (case units) of the quantity held fixed, for the "At …" select."""
    return q.grid_values(tables.summary, options.other)


def fixed_value(tables, options):
    values = fixed_choices(tables, options)
    if options.fixed is not None:
        return float(options.fixed)
    return values[0] if values else None


def flagged_x(points, tables, options):
    """x values (display units) of flagged points [(T_C, P_bar), ...] at the fixed value."""
    fixed = fixed_value(tables, options)
    out = []
    for T, P in points:
        x, other = (T, P) if options.x_axis == 'T' else (P, T)
        if fixed is not None and np.isclose(float(other), fixed, rtol=1e-9, atol=1e-9):
            out.append(options.convert(options.axis, x))
    return out


def map_items(tables, amount='mol', basis='dry'):
    """Choices of the T–P map: {item key: label}."""
    items = {'gas_carbon_percent': 'Carbon gasification efficiency [%]',
             'gibbs_J_per_mol_C': 'Gibbs energy [J / mol C]'}
    items.update({f'carbon:{p}': f'Carbon in {PHASE_TITLES[p].lower()} [% of feed carbon]'
                  for p in PANEL_ORDER})
    species = tables.species
    if species is not None and not species.empty:
        part = species[species.mol_per_mol_C.astype(float) > q.TINY]
        for phase, name in sorted(set(zip(part.phase, part.species))):
            items[f'species:{phase}:{name}'] = (f'{name}, {PHASE_TITLES.get(phase, phase).lower()} '
                                               f'[{q.amount_unit(amount, basis)}]')
    return items


def plot_values(kind, tables, options):
    """The long table a species plot draws from (T_C, P_bar, [phase,] species, value)."""
    if kind == 'total':
        return q.total_amounts(tables.species, options.amount, options.basis, options.s, options.a)
    if kind == 'phases':
        return q.amounts(tables.species, options.amount, options.basis, options.s, options.a)
    if kind == 'gas':
        return q.gas_composition(tables.species, options.composition, options.exclude_water)
    raise ValueError(f'{kind!r} has no species')


def largest_species(kind, tables, options, n=LARGEST):
    """The "Largest 8" button: the n largest species of the plot at the fixed value (no water)."""
    if kind not in SPECIES_PLOTS or tables.species is None or tables.species.empty:
        return []
    long = plot_values(kind, tables, options)
    fixed = fixed_value(tables, options)
    if fixed is not None:
        long = q.at(long, options.other, fixed)
    return q.largest(long, n, include_water=False)


def _selected(kind, tables, options):
    names = list(options.species) if options.species is not None else largest_species(kind, tables, options)
    names = [s for s in names if s != q.WATER]
    if options.include_water and not (kind == 'gas' and options.exclude_water):
        names.append(q.WATER)
    return names


# --- building -----------------------------------------------------------------------------------

def build_plot(kind, tables, options=None):
    """PlotSpec of one plot type (a key of PLOTS)."""
    options = options or PlotOptions()
    builders = {'total': amount_total, 'phases': amount_by_phase, 'gas': gas_composition,
                'cge': carbon_gasification_efficiency, 'carbon': carbon_in_phases, 'map': tp_map}
    return builders[kind](tables, options)


def _title(kind, unit, fixed_text=''):
    return ' · '.join(part for part in (PLOTS[kind], unit, fixed_text) if part)


def _no_points(kind):
    return PlotSpec(PLOTS[kind], [Panel('', empty_message=NO_POINTS)], kind, message=NO_POINTS)


def _ready(tables):
    return tables is not None and tables.summary is not None and not tables.summary.empty


def _colour(colours, name):
    if name not in colours:
        colours[name] = max(colours.values(), default=-1) + 1
    return colours[name]


def _series(table, options, colours, names=None, name_of=None):
    """Lines of a profile table (index = axis values in case units)."""
    x = [options.convert(options.axis, v) for v in table.index]
    out = []
    for name in (names if names is not None else table.columns):
        if name not in table:
            continue
        y = table[name].to_numpy(dtype=float)
        if options.log:
            y = np.where(y > LOG_FLOOR, y, np.nan)
        out.append(Series(x, y.tolist(), name_of(name) if name_of else name, _colour(colours, name), name))
    return out


def _cutoff(tables, options, x_values):
    """Dashed line at the aqueous cut-off on temperature axes within the plotted range."""
    if options.x_axis != 'T' or options.cutoff_C is None or not len(x_values):
        return []
    cutoff = float(options.cutoff_C)
    if not min(x_values) <= cutoff <= max(x_values):
        return []
    return [VLine(options.convert('T_C', cutoff), f"aqueous removed ≥ {options.text('T_C', cutoff)}")]


def _measured(tables, options, colours, names, phase=None):
    """Measured-data markers for the amount plots (same colour as the species)."""
    m = options.measured
    if m is None or m.table is None or m.table.empty:
        return []
    long = q.measured_amounts(m.table, m.quantity, tables.species, options.amount, options.basis,
                              options.s, options.a)
    if long.empty:
        return []
    fixed = fixed_value(tables, options)
    if options.x_axis == 'P':
        if 'P_bar' not in long or 'T_C' not in long:
            return []
        long = q.at(long, 'T_C', fixed)
    elif 'P_bar' in long and fixed is not None:
        long = q.at(long, 'P_bar', fixed)
    if phase is not None:
        long = long[long.phase == phase]
    out = []
    for (ph, name), part in long.groupby(['phase', 'species'], sort=False):
        if name not in names:
            continue
        part = part.sort_values(options.axis)
        y = part.value.to_numpy(dtype=float)
        if options.log:
            y = np.where(y > LOG_FLOOR, y, np.nan)
        label = f'{name} ({m.label})' if phase is not None else f'{name}, {PHASE_TITLES.get(ph, ph).lower()} ({m.label})'
        out.append(Series([options.convert(options.axis, v) for v in part[options.axis]], y.tolist(), label,
                          _colour(colours, name), f'measured:{ph}:{name}', None, 'square-open', False, True,
                          legend=False))
    if out:  # one legend entry for all measured points, in a neutral colour
        out.append(Series([], [], m.label, NEUTRAL, 'measured', None, 'square-open', False, True))
    return out


def _amount_setup(kind, tables, options):
    fixed = fixed_value(tables, options)
    colours = dict(options.colours) if options.colours is not None else run_colours(tables.species)
    unit = q.amount_unit(options.amount, options.basis)
    title = _title(kind, unit, options.text(options.other, fixed))
    return fixed, colours, unit, title


def amount_total(tables, options):
    """"Amount, all phases together": each species summed over the phases."""
    if not _ready(tables):
        return _no_points('total')
    fixed, colours, unit, title = _amount_setup('total', tables, options)
    names = _selected('total', tables, options)
    table = q.profile(plot_values('total', tables, options), tables.summary, options.axis, fixed, names)
    panel = Panel('', options.axis_label(options.axis), f'Amount [{unit}]', _series(table, options, colours),
                  options.log, key='total')
    panel.series += _measured(tables, options, colours, names)
    panel.vlines = _cutoff(tables, options, table.index)
    return PlotSpec(title, [panel], 'total')


def amount_by_phase(tables, options):
    """"Amount in each phase": 2 × 2 panels Gas, Aqueous, Solid, Organic liquid.

    A phase with no amount at any plotted point gets an empty panel with the message
    "No <phase> at any point".
    """
    if not _ready(tables):
        return _no_points('phases')
    fixed, colours, unit, title = _amount_setup('phases', tables, options)
    names = _selected('phases', tables, options)
    long = plot_values('phases', tables, options)
    rows = q.at(tables.summary, options.other, fixed) if fixed is not None else tables.summary
    points = list(zip(rows.T_C, rows.P_bar))
    panels = []
    for phase in PANEL_ORDER:
        title_p = PHASE_TITLES[phase]
        if q.phase_is_empty(tables.species, phase, points):
            panels.append(Panel(title_p, options.axis_label(options.axis), f'Amount [{unit}]', log=options.log,
                                empty_message=f'No {title_p.lower()} at any point', key=phase))
            continue
        part = long[long.phase == phase]
        keys = [s for s in names if s in set(part.species)]
        table = q.profile(part, tables.summary, options.axis, fixed, keys)
        panel = Panel(title_p, options.axis_label(options.axis), f'Amount [{unit}]',
                      _series(table, options, colours), options.log, key=phase)
        panel.series += _measured(tables, options, colours, names, phase)
        if phase == 'aqueous':
            panel.vlines = _cutoff(tables, options, table.index)
        panels.append(panel)
    return PlotSpec(title, panels, 'phases', rows=2, cols=2)


def gas_composition(tables, options):
    """Gas composition in mol % or mass %, optionally without water (dry gas)."""
    if not _ready(tables):
        return _no_points('gas')
    fixed = fixed_value(tables, options)
    colours = dict(options.colours) if options.colours is not None else run_colours(tables.species)
    unit = 'mol %' if options.composition == 'mol' else 'mass %'
    unit_title = f'{unit}, dry gas' if options.exclude_water else unit
    names = _selected('gas', tables, options)
    table = q.profile(plot_values('gas', tables, options), tables.summary, options.axis, fixed, names)
    panel = Panel('', options.axis_label(options.axis), f'Gas composition [{unit}]',
                  _series(table, options, colours), options.log, key='gas')
    return PlotSpec(_title('gas', unit_title, options.text(options.other, fixed)), [panel], 'gas')


def carbon_gasification_efficiency(tables, options):
    """Carbon gasification efficiency [%] = gas-phase carbon / feed carbon."""
    if not _ready(tables):
        return _no_points('cge')
    fixed = fixed_value(tables, options)
    long = q.carbon_gasification_efficiency(tables.summary).assign(species='cge')
    table = q.profile(long, tables.summary, options.axis, fixed, ['cge'])
    series = _series(table, options, {'cge': 0}, name_of=lambda _: PLOTS['cge'])
    panel = Panel('', options.axis_label(options.axis), f"{PLOTS['cge']} [%]", series, options.log, key='cge')
    return PlotSpec(_title('cge', '%', options.text(options.other, fixed)), [panel], 'cge')


def carbon_in_phases(tables, options):
    """Share of feed carbon in each phase [%] (phases that hold carbon)."""
    if not _ready(tables):
        return _no_points('carbon')
    fixed = fixed_value(tables, options)
    long = q.phase_carbon(tables.phases, tables.summary)
    table = q.profile(long, tables.summary, options.axis, fixed, list(PANEL_ORDER), key_column='phase')
    shown = [p for p in PANEL_ORDER if table[p].abs().max() > SHOWN]
    series = _series(table, options, dict(PHASE_COLOURS), shown, name_of=lambda p: PHASE_TITLES[p])
    panel = Panel('', options.axis_label(options.axis), 'Carbon [% of feed carbon]', series, options.log,
                  key='carbon')
    panel.vlines = _cutoff(tables, options, table.index)
    return PlotSpec(_title('carbon', '% of feed carbon', options.text(options.other, fixed)), [panel], 'carbon')


def tp_map(tables, options):
    """Colour map of one item over the temperature-pressure grid; failed points are empty."""
    if not _ready(tables):
        return _no_points('map')
    summary, item = tables.summary, options.map_item
    temperatures, pressures = q.grid_values(summary, 'T_C'), q.grid_values(summary, 'P_bar')
    label = map_items(tables, options.amount, options.basis).get(item, item)
    name, _, unit = label.rpartition(' [')
    unit = unit.rstrip(']')
    if item.startswith('species:'):
        phase, species = item[len('species:'):].split(':', 1)
        long = q.amounts(tables.species, options.amount, options.basis, options.s, options.a)
        long = long[(long.phase == phase) & (long.species == species)]
        fill = 0.0
    elif item.startswith('carbon:'):
        long = q.phase_carbon(tables.phases, summary)
        long = long[long.phase == item[len('carbon:'):]]
        fill = 0.0
    else:
        part = q.solved_rows(summary)
        long = part[['T_C', 'P_bar']].assign(value=part[item].astype(float))
        fill = None
    grid = long.pivot_table(index='P_bar', columns='T_C', values='value', aggfunc='sum') if len(long) else \
        pd.DataFrame(dtype=float)
    grid.index = grid.index.astype(float)
    grid.columns = grid.columns.astype(float)
    grid = grid.reindex(index=pressures, columns=temperatures).astype(float)
    solved = q.solved_rows(summary)
    done = {(float(T), float(P)) for T, P in zip(solved.T_C, solved.P_bar)}
    for P in pressures:
        for T in temperatures:
            if (T, P) not in done:
                grid.loc[P, T] = np.nan
            elif fill is not None and np.isnan(grid.loc[P, T]):
                grid.loc[P, T] = fill
    heatmap = Heatmap([options.number_format(options.convert('T_C', T)) for T in temperatures],
                      [options.number_format(options.convert('P_bar', P)) for P in pressures],
                      grid.to_numpy().tolist(), unit)
    panel = Panel('', options.axis_label('T_C'), options.axis_label('P_bar'), heatmap=heatmap, key='map')
    return PlotSpec(_title('map', name), [panel], 'map')
