"""Compare tab: one quantity from several runs in the same plot."""
import pandas as pd
from nicegui import ui

from .. import plotbuild as pb
from .. import quantities as q
from ..plotrender import render_plotly
from .common import OptionRadio, data_table, show_table
from .differences import setting_differences
from .results_tab import (
    PRESSURE_UNITS,
    TEMPERATURE_UNITS,
    model_molar_masses,
    run_feed_fractions,
    run_measured,
    unit_options,
)
from .save_plot import SavePlotDialog

CAPTION = 'One quantity from several runs in the same plot.'
NO_MEASURED = 'None'
SAME_SETTINGS = 'The runs have the same settings.'


class CompareTab:
    def __init__(self, state):
        self.state = state
        self.settings = state.settings
        self.loaded = []  # CompareRun per selected run
        self.spec = None
        self.save_dialog = None
        self._masses = None
        self._quiet = False
        self.t_unit = self.settings.get('temperature_unit')
        self.p_unit = self.settings.get('pressure_unit')
        ui.label(CAPTION).classes('text-caption text-grey-7')
        with ui.row().classes('w-full items-end gap-3'):
            self.runs = ui.select(state.sources(), label='Runs to compare', multiple=True, with_input=True,
                                  on_change=self.load).classes('w-[40rem]').props('use-chips')
            ui.button(icon='refresh', on_click=self.refresh_sources).props('flat round').tooltip('Refresh the list')
        with ui.row().classes('w-full items-end gap-4'):
            self.key = ui.select({}, label='Quantity', with_input=True, on_change=self.draw).classes('w-80')
            with ui.column().classes('gap-0'):
                ui.label('Amount as').classes('text-caption')
                self.amount = ui.toggle(q.AMOUNTS, value='mol', on_change=self.draw).props('dense no-caps')
            with ui.column().classes('gap-0'):
                ui.label('per kg of').classes('text-caption')
                self.basis = OptionRadio(q.BASIS_OPTIONS, value='dry', on_change=self.draw).props('dense inline')
            with ui.column().classes('gap-0'):
                ui.label('Temperature unit').classes('text-caption')
                self.t_toggle = ui.toggle(TEMPERATURE_UNITS, value=self.t_unit,
                                          on_change=self._unit_changed).props('dense no-caps')
            with ui.column().classes('gap-0'):
                ui.label('Pressure unit').classes('text-caption')
                self.p_toggle = ui.toggle(PRESSURE_UNITS, value=self.p_unit,
                                          on_change=self._unit_changed).props('dense no-caps')
            self.log = ui.checkbox('Log scale', on_change=self.draw).props('dense')
            self.measured = ui.select({'': NO_MEASURED}, label='Measured data', value='',
                                      on_change=self.draw).classes('w-72')
        with ui.card().classes('w-full'):
            with ui.row().classes('w-full items-center'):
                self.message = ui.label().classes('text-grey-7')
                ui.space()
                self.save_button = ui.button('Save plot…', icon='save', on_click=self.open_save) \
                    .props('outline no-caps')
            self.plot = ui.plotly({}).classes('w-full h-[32rem]')
        ui.label('Settings that differ').classes('text-h6')
        self.same = ui.label(SAME_SETTINGS).classes('text-grey-7')
        self.differences = data_table(pd.DataFrame(), rows_per_page=20)
        self.draw()
        ui.timer(5.0, self.refresh_sources)

    def refresh_sources(self):
        sources = self.state.sources()
        if sources != self.runs.options:
            self.runs.set_options(sources, value=[v for v in self.runs.value or [] if v in sources])

    def molar_masses(self):
        if self._masses is None:
            model = getattr(self.state, 'model', None)
            self._masses = model_molar_masses(model) if model is not None else {}
        return self._masses

    def load(self, _=None):
        self.loaded = []
        for source in self.runs.value or []:
            try:
                tables = self.state.tables(source)
            except Exception as exc:
                ui.notify(f'Could not read {source}: {exc}', type='negative')
                continue
            if tables is None:
                continue
            try:
                fractions = run_feed_fractions(tables, self.molar_masses())
            except Exception:
                fractions = None
            s, a = fractions if fractions is not None else (None, 0.0)
            self.loaded.append(pb.CompareRun(tables.label, tables, s, a))
        self._quiet = True
        try:
            items = pb.compare_items(self.loaded)
            self.key.set_options(items, value=self.key.value if self.key.value in items else pb.CGE)
            measured = {'': NO_MEASURED}
            for run in self.loaded:
                found = run_measured(run.tables)
                if found is not None:
                    measured[run.label] = f'{found.label} (from {run.label})'
            self.measured.set_options(measured, value=self.measured.value if self.measured.value in measured else '')
            missing = [run.label for run in self.loaded if run.s is None]
            self.basis.set_disabled({'total': q.NO_TOTAL_MASS} if missing else {})
            if missing and self.basis.value == 'total':
                self.basis.set_value('dry')
        finally:
            self._quiet = False
        self.fill_differences()
        self.draw()

    def _unit_changed(self, _=None):
        if self._quiet:
            return
        self.t_unit = self.t_toggle.value or self.t_unit
        self.p_unit = self.p_toggle.value or self.p_unit
        self.fill_differences()
        self.draw()

    def fill_differences(self):
        table = setting_differences({run.label: run.tables.case for run in self.loaded}, self.t_unit, self.p_unit)
        show_table(self.differences, table)
        self.same.set_visibility(len(self.loaded) > 1 and table.empty)

    def options(self):
        measured_run = next((run for run in self.loaded if run.label == self.measured.value), None)
        measured = run_measured(measured_run.tables) if measured_run is not None else None
        return unit_options(self.t_unit, self.p_unit, amount=self.amount.value, basis=self.basis.value,
                            log=bool(self.log.value), measured=measured), measured_run

    def draw(self, _=None):
        if self._quiet:
            return
        key = self.key.value or pb.CGE
        for element in (self.amount, self.basis):
            element.set_enabled(key != pb.CGE)
        self.measured.set_enabled(key != pb.CGE)
        if not self.loaded:
            self.spec = None
            self.message.set_text('Choose two or more runs.')
            self.plot.update_figure({})
            self.save_button.disable()
            return
        options, measured_run = self.options()
        try:
            spec = pb.compare_runs(self.loaded, key, options, measured_run)
            figure = render_plotly(spec, self.settings.get('color_theme'), self.settings.get('line_width'))
        except Exception as exc:
            self.spec = None
            self.message.set_text(f'Cannot draw this plot: {exc}')
            self.plot.update_figure({})
            self.save_button.disable()
            return
        self.spec = spec
        self.message.set_text(f'{len(self.loaded)} run{"s" if len(self.loaded) != 1 else ""}')
        self.plot.update_figure(figure)
        self.save_button.set_enabled(not spec.message)

    def open_save(self):
        if self.spec is None:
            ui.notify('Choose two or more runs.')
            return
        if self.save_dialog is not None:
            self.save_dialog.dialog.delete()
        stem = 'compare_' + (self.key.value or pb.CGE).replace(':', '_').replace(' ', '_')
        self.save_dialog = SavePlotDialog(self.settings, self.spec, self.settings.get('color_theme'),
                                          self.settings.get('line_width'), stem, layout='one').open()
