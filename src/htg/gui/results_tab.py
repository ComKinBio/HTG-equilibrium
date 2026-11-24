"""Results tab: plots and tables of one run.

Left: the settings panel. Right: the banner area, the plot card and the tables
card. Every plot is built as a PlotSpec (htg.plotbuild) and drawn with Plotly;
"Save plot…" draws the same spec with Matplotlib.

Points with more than one result (WP6): ``show_banners()`` fills the banner
above the plot; the dialog (htg.gui.multiple_results) chooses the result each
point uses (choices.toml, applied by RunTables.with_choices) and which results
comparison mode draws; ``decorate(spec, options)`` turns a spec into the
comparison spec (plotspec.compare_results).
"""
import io
import tomllib
import zipfile
from dataclasses import replace
from pathlib import Path

import pandas as pd
from nicegui import ui

from .. import plotbuild as pb
from .. import quantities as q
from .. import units
from ..case import Feed
from ..feed_v2 import feed_fractions, feed_v2_to_model, is_v2_feed, water_molar_mass
from ..plotrender import render_plotly
from ..plotspec import (
    LINE_WIDTHS,
    THEME_LABELS,
    Series,
    compare_results,
    palette,
    resolve_style,
)
from ..results import write_choices
from .common import (
    OptionRadio,
    data_table,
    palette_strip,
    save_bytes,
    show_table,
    swatch,
)
from .multiple_results import (
    MultipleResultsDialog,
    banner_text,
    comparing_text,
    results_at,
)
from .save_plot import SavePlotDialog

TEMPERATURE_UNITS = {key: units.LABELS[key] for key in units.TEMPERATURE_UNITS}
PRESSURE_UNITS = {key: units.LABELS[key] for key in units.PRESSURE_UNITS}
X_AXES = {'T': 'Temperature', 'P': 'Pressure'}
COMPOSITIONS = {'mol': 'mol %', 'mass': 'mass %'}
MAP_NEEDS_GRID = 'The T–P map needs at least 2 temperatures and 2 pressures.'
TABLES = ('summary', 'species', 'phases')


# ---------------------------------------------------------------- run information
def run_feed(tables):
    """The model Feed of a run: from the Case it was run with, else from its case.toml."""
    case = getattr(tables, 'case', None)
    if case is not None:
        return case.feed
    path = getattr(tables, 'path', None)
    if path is None or not (Path(path) / 'case.toml').exists():
        return None
    try:
        with open(Path(path) / 'case.toml', 'rb') as f:
            feed = tomllib.load(f).get('feed')
        if is_v2_feed(feed):
            return feed_v2_to_model(feed)
        return Feed(**feed) if isinstance(feed, dict) else None
    except Exception:  # an unreadable or invalid file: no feed information
        return None


def run_feed_fractions(tables, molar_masses=None):
    """(s, a) of a run's feed (htg.feed_v2.feed_fractions), or None without a total mass.

    s = solids / total feed, a = ash / solids. Old-format and version-2 cases both
    work, through the model Feed. None disables the "total feed" basis.
    molar_masses  {species: g/mol} for feeds given as species amounts
    """
    feed = run_feed(tables) if tables is not None else None
    if feed is None:
        return None
    masses = dict(molar_masses or {})
    species = getattr(tables, 'species', None)
    if species is not None and not species.empty:
        masses = {**q.molar_masses(species), **masses}
    masses.setdefault(q.WATER, water_molar_mass())
    try:
        s, a = feed_fractions(feed, masses.__getitem__)
    except Exception:
        return None
    return None if s is None else (s, a)


def run_measured(tables):
    """Measured data saved with a run (its comparison file), or None."""
    case = getattr(tables, 'case', None)
    settings = getattr(case, 'comparison', None) if case is not None else None
    if tables is None or tables.comparison is None or not settings:
        return None
    return pb.Measured(tables.comparison, settings.get('quantity', 'mol_per_kg_dry_feed'),
                       settings.get('label') or 'Measured')


def run_cutoff(tables):
    case = getattr(tables, 'case', None)
    return getattr(case, 'aqueous_cutoff_C', None) if case is not None else None


def unit_options(t_unit, p_unit, **changes):
    """PlotOptions with display units (keys of htg.units, e.g. 'K', 'MPa')."""
    return pb.PlotOptions(t_unit=units.label(t_unit), p_unit=units.label(p_unit),
                          t_convert=lambda v: units.from_file(v, t_unit),
                          p_convert=lambda v: units.from_file(v, p_unit),
                          t_format=lambda v: units.format_number(v, t_unit),
                          p_format=lambda v: units.format_number(v, p_unit), **changes)


def model_molar_masses(model):
    """{species: g/mol} of the model's catalogue (for feeds given in mol)."""
    from ..formulas import molar_mass, parse_formula

    out = {}
    try:
        catalogue = model.species()
    except Exception:
        return out
    for name, formula in zip(catalogue.Species, catalogue.Formula):
        try:
            out[str(name)] = molar_mass(parse_formula(str(formula)))
        except Exception:
            continue
    return out


def heading(text):
    return ui.label(text).classes('text-weight-medium q-mt-sm')


# ---------------------------------------------------------------- the tab
class ResultsTab:
    def __init__(self, state):
        self.state = state
        self.settings = state.settings
        self.tables = None
        self.fractions = None
        self.colours = {}
        self.selected = None  # species ticked by the user; None = the largest 8 of the plot
        self.spec = None
        self.shown = None  # (source, done points) drawn last
        self.save_dialog = None
        self._masses = None
        self._quiet = False
        self._units_changed = False
        self.t_unit = self.settings.get('temperature_unit')
        self.p_unit = self.settings.get('pressure_unit')
        self.species_checks = {}
        self._species_names = None
        self.comparing = False  # comparison mode (WP6)
        self.compare_show = {}  # {(T, P): {result: shown when comparing}}
        self.results_dialog = None

        with ui.row().classes('w-full no-wrap items-start gap-4'):
            with ui.card().classes('w-80 shrink-0 gap-1'):
                self._settings_panel()
            with ui.column().classes('grow min-w-0 gap-3'):
                self.banner = ui.column().classes('w-full gap-2')
                self._plot_card()
                self._tables_card()
        self.dialog_host = ui.element('div')  # dialogs live here, not in the banner that is rebuilt
        self.controls()
        self.draw()
        ui.timer(2.0, self.follow)
        ui.timer(5.0, self.refresh_sources)  # new runs and results folders

    # ------------------------------------------------------------ layout
    def _settings_panel(self):
        with ui.row().classes('w-full items-end no-wrap gap-1'):
            self.source = ui.select(self.state.sources(), label='Run', with_input=True,
                                    on_change=self.load).classes('grow')
            ui.button(icon='refresh', on_click=self.refresh_sources).props('flat round dense') \
                .tooltip('Refresh the list')
        self.status = ui.label().classes('text-caption text-grey-7')
        heading('Plot')
        self.kind = OptionRadio(pb.PLOTS, value='phases', on_change=self.draw).props('dense')
        self.map_item = ui.select({}, label='Map quantity', with_input=True, on_change=self.draw).classes('w-full')
        with ui.column().classes('w-full gap-1') as self.amount_box:
            heading('Amount as')
            self.amount = ui.toggle(q.AMOUNTS, value='mol', on_change=self.draw).props('dense no-caps')
            ui.label('per kg of').classes('text-caption')
            self.basis = OptionRadio(q.BASIS_OPTIONS, value='dry', on_change=self.draw).props('dense')
        with ui.column().classes('w-full gap-1') as self.gas_box:
            heading('Gas composition only')
            self.composition = ui.toggle(COMPOSITIONS, value='mol', on_change=self.draw).props('dense no-caps')
            self.exclude_water = ui.checkbox('Exclude water (dry gas)', on_change=self.draw).props('dense')
        heading('x axis')
        with ui.row().classes('w-full items-center gap-2'):
            self.x_axis = ui.toggle(X_AXES, value='T', on_change=self._axis_changed).props('dense no-caps')
            self.x_unit = ui.toggle(TEMPERATURE_UNITS, value=self.t_unit,
                                    on_change=self._x_unit_changed).props('dense no-caps')
        with ui.row().classes('w-full items-end no-wrap gap-2'):
            self.fixed = ui.select({}, label='At pressure', on_change=self.draw).classes('grow')
            self.fixed_unit = ui.toggle(PRESSURE_UNITS, value=self.p_unit,
                                        on_change=self._fixed_unit_changed).props('dense no-caps')
        with ui.column().classes('w-full gap-1') as self.species_box:
            heading('Species')
            with ui.row().classes('w-full items-end no-wrap gap-2'):
                ui.button(f'Largest {pb.LARGEST}', on_click=self.largest).props('outline dense no-caps')
                self.filter = ui.input('Filter', on_change=self.fill_species).props('dense clearable') \
                    .classes('grow')
            self.water = ui.checkbox('water', on_change=self.draw).props('dense')
            self.species_list = ui.column().classes('w-full gap-0 max-h-72 overflow-auto')
        self.log = ui.checkbox('Log scale', on_change=self.draw).props('dense')

    def _plot_card(self):
        with ui.card().classes('w-full'):
            with ui.row().classes('w-full items-center gap-2'):
                self.message = ui.label().classes('text-grey-7')
                ui.space()
                theme = self.settings.get('color_theme')
                self.theme = ui.select(THEME_LABELS, label='Colours',
                                       value=theme if theme in THEME_LABELS else 'Default',
                                       on_change=self._theme_changed).classes('w-40')
                self.strip = ui.row().classes('gap-0 no-wrap items-center')
                width = self.settings.get('line_width')
                self.line = ui.select({w: f'{w:g} px' for w in LINE_WIDTHS}, label='Line',
                                      value=width if width in LINE_WIDTHS else 2,
                                      on_change=self._line_changed).classes('w-24')
                self.save_button = ui.button('Save plot…', icon='save', on_click=self.open_save).props('outline no-caps')
            palette_strip(self.strip, palette(self.theme.value))
            self.plot = ui.plotly({}).classes('w-full h-[40rem]')

    def _tables_card(self):
        with ui.card().classes('w-full'):
            with ui.row().classes('w-full items-center'):
                ui.label('Tables').classes('text-subtitle1 text-weight-medium')
                ui.space()
                ui.button('Download all as CSV', icon='download', on_click=self.download).props('flat dense no-caps') \
                    .tooltip('One ZIP file with summary.csv, species.csv and phases.csv')
            with ui.tabs().classes('w-full') as tabs:
                summary_tab = ui.tab('Summary')
                species_tab = ui.tab('Species')
                phases_tab = ui.tab('Phases')
            with ui.tab_panels(tabs, value=summary_tab).classes('w-full'):
                with ui.tab_panel(summary_tab):
                    self.summary_table = data_table(pd.DataFrame())
                with ui.tab_panel(species_tab):
                    self.species_table = data_table(pd.DataFrame())
                with ui.tab_panel(phases_tab):
                    self.phases_table = data_table(pd.DataFrame())

    # ------------------------------------------------------------ runs
    def refresh_sources(self):
        sources = self.state.sources()
        if sources != self.source.options:
            self.source.set_options(sources, value=self.source.value if self.source.value in sources else None)

    def select(self, source):
        self.refresh_sources()
        if self.source.value == source:
            self.load()
        else:
            self.source.set_value(source)

    def _done(self, source):
        run = self.state.run(source)
        return run.done_points if run else None

    def follow(self):
        """Redraw while a run of this session is still producing points."""
        source = self.source.value
        run = self.state.run(source)
        if run is not None and (source, run.done_points) != self.shown:
            self.refresh_sources()
            self.load()

    def molar_masses(self):
        if self._masses is None:
            model = getattr(self.state, 'model', None)
            self._masses = model_molar_masses(model) if model is not None else {}
        return self._masses

    def load(self, _=None):
        source = self.source.value
        same_run = self.shown is not None and self.shown[0] == source
        try:
            self.tables = self.state.tables(source)
        except Exception as exc:
            self.tables = None
            ui.notify(f'Could not read the results: {exc}', type='negative', multi_line=True)
        tables = self.tables
        if not same_run:
            self.selected = None
            self.comparing = False
            self.compare_show = {}
        if not self._units_changed:
            self.t_unit = self.settings.get('temperature_unit')
            self.p_unit = self.settings.get('pressure_unit')
        self.fractions = None
        self.colours = {}
        self._quiet = True
        try:
            if tables is not None:
                try:
                    self.fractions = run_feed_fractions(tables, self.molar_masses())
                except Exception:
                    self.fractions = None
                self.colours = pb.run_colours(tables.species)
                if not same_run and len(tables.temperatures) == 1 and len(tables.pressures) > 1:
                    self.x_axis.set_value('P')
            self.basis.set_disabled({} if self.fractions is not None else {'total': q.NO_TOTAL_MASS})
            if self.basis.value == 'total' and self.fractions is None:
                self.basis.set_value('dry')
            grid = tables is not None and pb.map_available(tables.summary)
            self.kind.set_disabled({} if grid else {'map': MAP_NEEDS_GRID})
            if self.kind.value == 'map' and not grid:
                self.kind.set_value('phases')
            self._unit_toggles()
            self._fixed_options()
            self._map_options()
        finally:
            self._quiet = False
        self.shown = (source, self._done(source))
        self._species_names = None
        self.show_banners()
        self.draw()
        self.fill_tables()

    # ------------------------------------------------------------ points with more than one result
    def view_tables(self):
        """The tables shown in plots and tables (with the chosen results, see choices.toml)."""
        return self.tables

    def compared(self):
        """({result: points drawn for it}, results in the line-style order) of comparison mode."""
        tables = self.tables
        per_label = {}
        for point in tables.flagged_points if tables is not None else []:
            shown = self.compare_show.get(point, {})
            for label in results_at(tables.alternatives, point):
                if label != 'A' and shown.get(label, True):
                    per_label.setdefault(label, []).append(point)
        return per_label, ['A', *sorted(per_label)]

    def decorate(self, spec, options):
        """The comparison spec in comparison mode: result A solid, the other ticked results
        dashed, dotted, … around the points that have them (plotspec.compare_results)."""
        if not self.comparing or spec.message or self.kind.value == 'map':
            return spec
        per_label, _ = self.compared()
        if not per_label:
            return spec
        kind = self.kind.value
        options = replace(options, species=self.shown_species(options))
        base = pb.build_plot(kind, self.tables.with_choices({}), options)
        others = {label: pb.build_plot(kind, self.tables.with_choices({p: label for p in points}), options)
                  for label, points in per_label.items()}
        flagged = {label: pb.flagged_x(points, self.tables, options) for label, points in per_label.items()}
        return compare_results(base, others, flagged, descriptions={'A': 'lowest Gibbs energy'})

    def show_banners(self):
        """The banner above the plot: points with more than one result, or comparison mode."""
        self.banner.clear()
        tables = self.tables
        flagged = tables.flagged_points if tables is not None else []
        if not flagged:
            self.comparing = False
            return
        with self.banner, ui.card().classes('w-full bg-amber-1 q-py-sm').props('flat bordered'):
            with ui.row().classes('w-full items-center no-wrap gap-2'):
                ui.icon('info', color='amber-9')
                if self.comparing:
                    per_label, labels = self.compared()
                    points = {p for ps in per_label.values() for p in ps}
                    ui.label(comparing_text(len(points), labels)).classes('grow')
                    ui.button('Change', on_click=self.open_results).props('flat dense no-caps')
                    ui.button('Stop comparing', on_click=self.stop_comparing).props('flat dense no-caps')
                else:
                    ui.label(banner_text(len(flagged), len(tables.summary))).classes('grow')
                    ui.button('Compare them', on_click=self.open_results).props('flat dense no-caps')

    def open_results(self):
        if self.tables is None or not self.tables.flagged_points:
            return
        if self.results_dialog is not None and not self.results_dialog.dialog.is_deleted:
            self.results_dialog.dialog.delete()
        total = self.fractions is not None
        with self.dialog_host:
            self.results_dialog = MultipleResultsDialog(
                self.tables, self.compare_show, self.apply_choices, self.start_comparing, self.t_unit,
                self.p_unit, factor=self.fractions[0] if total else 1.0,
                unit_text='mol / kg total feed' if total else 'mol / kg dry feed').open()

    def apply_choices(self, choices):
        """Use the chosen results in plots and tables and save them in choices.toml."""
        folder = self.tables.path
        if folder is not None:
            Path(folder).mkdir(parents=True, exist_ok=True)
            write_choices(folder, choices)
        else:
            ui.notify('This run has no results folder, so the choices are not saved.')
        self.tables = self.tables.with_choices(choices)
        self._species_names = None
        self.show_banners()
        self.draw()
        self.fill_tables()

    def start_comparing(self):
        if not self.compared()[0]:
            ui.notify('Tick at least one result other than A to compare.')
            return
        self.comparing = True
        self.show_banners()
        self.draw()

    def stop_comparing(self):
        self.comparing = False
        self.show_banners()
        self.draw()

    # ------------------------------------------------------------ units and choices
    def _unit_toggles(self):
        """Unit switches for the x axis and the fixed value, from t_unit and p_unit."""
        temperature_x = self.x_axis.value == 'T'
        x_units, x_value = (TEMPERATURE_UNITS, self.t_unit) if temperature_x else (PRESSURE_UNITS, self.p_unit)
        f_units, f_value = (PRESSURE_UNITS, self.p_unit) if temperature_x else (TEMPERATURE_UNITS, self.t_unit)
        self.x_unit.set_options(x_units, value=x_value)
        self.fixed_unit.set_options(f_units, value=f_value)
        self.fixed.label = 'At pressure' if temperature_x else 'At temperature'

    def _fixed_options(self):
        """The run's values of the quantity held fixed, in the current unit."""
        tables = self.view_tables()
        if tables is None or tables.summary.empty:
            self.fixed.set_options({}, value=None)
            return
        column = 'P_bar' if self.x_axis.value == 'T' else 'T_C'
        unit = self.p_unit if column == 'P_bar' else self.t_unit
        values = q.grid_values(tables.summary, column)
        options = {v: f'{units.display(v, unit)} {units.label(unit)}' for v in values}
        value = self.fixed.value if self.fixed.value in options else (values[0] if values else None)
        self.fixed.set_options(options, value=value)

    def _map_options(self):
        tables = self.view_tables()
        items = pb.map_items(tables, self.amount.value, self.basis.value) if tables is not None else {}
        value = self.map_item.value if self.map_item.value in items else 'gas_carbon_percent'
        self.map_item.set_options(items, value=value if items else None)

    def _axis_changed(self, _=None):
        if self._quiet:
            return
        self._quiet = True
        try:
            self._unit_toggles()
            self.fixed.set_value(None)
            self._fixed_options()
        finally:
            self._quiet = False
        self.draw()

    def _set_unit(self, unit):
        if unit is None:
            return
        if unit in units.TEMPERATURE_UNITS:
            self.t_unit = unit
        else:
            self.p_unit = unit
        self._units_changed = True
        self._quiet = True
        try:
            self._fixed_options()
        finally:
            self._quiet = False
        self.draw()

    def _x_unit_changed(self, e):
        if not self._quiet:
            self._set_unit(e.value)

    def _fixed_unit_changed(self, e):
        if not self._quiet:
            self._set_unit(e.value)

    def _theme_changed(self, e):
        if e.value in THEME_LABELS:
            self.settings.set('color_theme', e.value)
        palette_strip(self.strip, palette(self.theme.value))
        self._species_names = None
        self.draw()

    def _line_changed(self, e):
        if e.value in LINE_WIDTHS:
            self.settings.set('line_width', e.value)
        self.draw()

    def options(self):
        """PlotOptions from the settings panel."""
        tables = self.view_tables()
        s, a = self.fractions if self.fractions is not None else (None, 0.0)
        basis = self.basis.value if not (self.basis.value == 'total' and s is None) else 'dry'
        return unit_options(
            self.t_unit, self.p_unit, x_axis=self.x_axis.value, fixed=self.fixed.value, amount=self.amount.value,
            basis=basis, s=s, a=a, species=self.selected, include_water=bool(self.water.value),
            log=bool(self.log.value), composition=self.composition.value,
            exclude_water=bool(self.exclude_water.value), map_item=self.map_item.value or 'gas_carbon_percent',
            cutoff_C=run_cutoff(tables), measured=run_measured(tables), colours=self.colours or None)

    def controls(self):
        kind = self.kind.value
        self.map_item.set_visibility(kind == 'map')
        species_map = kind == 'map' and (self.map_item.value or '').startswith('species:')
        for element in (self.amount, self.basis):
            element.set_enabled(kind in pb.AMOUNT_PLOTS or species_map)
        for element in (self.composition, self.exclude_water):
            element.set_enabled(kind == 'gas')
        self.x_axis.set_enabled(kind != 'map')
        self.fixed.set_enabled(kind != 'map')
        species = kind in pb.SPECIES_PLOTS
        self.species_box.set_visibility(species)
        if self.tables is not None:
            solved = len(q.solved_rows(self.tables.summary))
            self.status.set_text(f'{solved} of {len(self.tables.summary)} points solved')
        else:
            self.status.set_text('')

    # ------------------------------------------------------------ species list
    def shown_species(self, options=None):
        """Species drawn now (without water)."""
        if self.selected is not None:
            return [s for s in self.selected if s != q.WATER]
        tables = self.view_tables()
        if tables is None:
            return []
        return pb.largest_species(self.kind.value, tables, options or self.options())

    def fill_species(self, _=None):
        """Rebuild the species checkboxes (run, filter or colour theme changed)."""
        tables = self.view_tables()
        names = [] if tables is None or tables.species.empty else \
            [n for n in pb.species_choices(tables.species) if n != q.WATER]
        text = (self.filter.value or '').strip().casefold()
        if text:
            names = [n for n in names if text in n.casefold()]
        shown = set(self.shown_species())
        self.species_list.clear()
        self.species_checks = {}
        theme = self.theme.value
        with self.species_list:
            for name in names:
                colour, dash = resolve_style(Series([], [], name, colour=self.colours.get(name, 0)), theme)
                with ui.row().classes('items-center no-wrap gap-1'):
                    check = ui.checkbox(name, value=name in shown,
                                        on_change=lambda e, n=name: self.tick(n, e.value)).props('dense')
                    swatch(colour, dash)
                self.species_checks[name] = check
            if not names:
                ui.label('No species match.' if text else 'No species.').classes('text-caption text-grey-7')
        self._species_names = (tuple(names), theme, text)

    def sync_species(self, options=None):
        if self._species_names is None:
            self.fill_species()
            return
        shown = set(self.shown_species(options))
        self._quiet = True
        try:
            for name, check in self.species_checks.items():
                if check.value != (name in shown):
                    check.set_value(name in shown)
        finally:
            self._quiet = False

    def tick(self, name, value):
        if self._quiet:
            return
        current = self.shown_species()
        if value and name not in current:
            current.append(name)
        elif not value:
            current = [s for s in current if s != name]
        self.selected = current
        self.draw()

    def largest(self):
        self.selected = None
        self.draw()

    # ------------------------------------------------------------ drawing
    def draw(self, _=None):
        if self._quiet:
            return
        self.controls()
        tables = self.view_tables()
        if tables is None:
            self.spec = None
            self.message.set_text('Choose a run.')
            self.plot.update_figure({})
            self.save_button.disable()
            self.sync_species()
            return
        self._quiet = True
        try:
            self._map_options()
        finally:
            self._quiet = False
        options = self.options()
        try:
            spec = self.decorate(pb.build_plot(self.kind.value, tables, options), options)
            figure = render_plotly(spec, self.theme.value, self.line.value)
        except Exception as exc:
            self.spec = None
            self.message.set_text(f'Cannot draw this plot: {exc}')
            self.plot.update_figure({})
            self.save_button.disable()
            return
        self.spec = spec
        self.message.set_text(spec.message)
        self.plot.classes(replace=f'w-full {"h-[44rem]" if spec.rows > 1 else "h-[32rem]"}')
        self.plot.update_figure(figure)
        self.save_button.set_enabled(not spec.message)
        self.sync_species(options)

    def stem(self):
        tables = self.tables
        name = getattr(tables.case, 'name', None) if tables is not None and tables.case is not None else None
        return name or 'results'

    def open_save(self):
        if self.spec is None:
            ui.notify('Choose a run first.')
            return
        if self.save_dialog is not None and not self.save_dialog.dialog.is_deleted:
            self.save_dialog.dialog.delete()
        with self.dialog_host:
            self.save_dialog = SavePlotDialog(self.settings, self.spec, self.theme.value, self.line.value,
                                              f'{self.stem()}_{self.kind.value}').open()

    # ------------------------------------------------------------ tables
    def fill_tables(self, _=None):
        tables = self.view_tables()
        empty = pd.DataFrame()
        show_table(self.summary_table, tables.summary if tables else empty)
        show_table(self.species_table, tables.species if tables else empty)
        show_table(self.phases_table, tables.phases if tables else empty)

    def tables_zip(self):
        """summary.csv, species.csv and phases.csv in one ZIP file (bytes)."""
        tables = self.view_tables()
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, 'w', zipfile.ZIP_DEFLATED) as archive:
            for name in TABLES:
                archive.writestr(f'{name}.csv', getattr(tables, name).to_csv(index=False))
        return buffer.getvalue()

    async def download(self):
        if self.view_tables() is None:
            ui.notify('Choose a run first.')
            return
        where = await save_bytes(self.tables_zip(), f'{self.stem()}_tables.zip')
        if where is not None and Path(str(where)).is_absolute():
            ui.notify(f'Saved {where}')
