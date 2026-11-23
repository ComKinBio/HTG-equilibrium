"""Results tab: plots and tables of one run."""
import pandas as pd
from nicegui import ui

from . import plots
from .common import data_table, show_table


VIEWS = {'temperature': 'Species vs temperature', 'pressure': 'Species vs pressure',
         'map': 'T–P map', 'carbon': 'Carbon allocation'}


class ResultsTab:
    def __init__(self, state):
        self.state = state
        self.tables = None
        self.shown = None  # (source, done_points) drawn last
        with ui.row().classes('w-full items-end gap-3'):
            self.source = ui.select(state.sources(), label='Run', with_input=True,
                                    on_change=self.load).classes('w-[28rem]')
            ui.button(icon='refresh', on_click=self.refresh_sources).props('flat round').tooltip('Refresh the list')
            self.view = ui.toggle(VIEWS, value='temperature', on_change=self.draw)
        with ui.row().classes('w-full items-end gap-3'):
            self.group = ui.select(plots.GROUPS, label='Group', value='gas', on_change=self.draw).classes('w-44')
            self.quantity = ui.select(plots.QUANTITIES, label='Quantity', value='mol_per_kg_dry_feed',
                                      on_change=self.draw).classes('w-56')
            self.top = ui.number('Largest N species', value=12, min=1, step=1, precision=0,
                                 on_change=self.draw).classes('w-36')
            self.log = ui.switch('Log scale', value=True, on_change=self.draw)
            self.pressures = ui.select([], label='Pressures [bar]', multiple=True,
                                       on_change=self.draw).classes('w-56').props('use-chips')
            self.temperature = ui.select([], label='Temperature [°C]', on_change=self.draw).classes('w-40')
            self.item = ui.select({}, label='Map quantity', with_input=True, on_change=self.draw).classes('w-72')
        self.message = ui.label().classes('text-grey-7')
        self.plot = ui.plotly(plots.go.Figure()).classes('w-full h-[34rem]')
        with ui.row().classes('items-center gap-2'):
            ui.label('Download CSV:')
            for name in ('summary', 'species', 'phases'):
                ui.button(name, icon='download', on_click=lambda n=name: self.download(n)).props('flat dense')
        with ui.tabs().classes('w-full') as tabs:
            summary_tab = ui.tab('Summary')
            species_tab = ui.tab('Species')
            phases_tab = ui.tab('Phases')
        with ui.tab_panels(tabs, value=summary_tab).classes('w-full'):
            with ui.tab_panel(summary_tab):
                self.summary_table = data_table(pd.DataFrame())
            with ui.tab_panel(species_tab):
                self.species_filter = ui.input('Filter species', on_change=self.fill_tables).classes('w-64')
                self.species_table = data_table(pd.DataFrame())
            with ui.tab_panel(phases_tab):
                self.phases_table = data_table(pd.DataFrame())
        self.controls()
        ui.timer(2.0, self.follow)
        ui.timer(5.0, self.refresh_sources)  # new runs and results folders

    def refresh_sources(self):
        sources = self.state.sources()
        if sources != self.source.options:
            self.source.set_options(sources, value=self.source.value if self.source.value in sources else None)

    def select(self, source):
        self.refresh_sources()
        self.source.set_value(source)

    def load(self, _=None):
        source = self.source.value
        try:
            self.tables = self.state.tables(source)
        except Exception as exc:
            self.tables = None
            ui.notify(f'Could not read the results: {exc}', type='negative', multi_line=True)
        tables = self.tables
        if tables is not None:
            pressures = [float(p) for p in tables.pressures]
            temperatures = [float(t) for t in tables.temperatures]
            keep = [p for p in (self.pressures.value or []) if p in pressures]
            self.pressures.set_options({p: f'{p:g}' for p in pressures}, value=keep or pressures[:6])
            self.temperature.set_options({t: f'{t:g}' for t in temperatures},
                                         value=self.temperature.value if self.temperature.value in temperatures
                                         else (temperatures[-1] if temperatures else None))
            items = plots.map_items(tables)
            self.item.set_options(items, value=self.item.value if self.item.value in items else 'gas_carbon_percent')
            if self.view.value == 'temperature' and len(temperatures) == 1 and len(pressures) > 1:
                self.view.set_value('pressure')
        self.shown = (source, self._done(source))
        self.draw()
        self.fill_tables()

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

    def controls(self):
        view = self.view.value
        self.group.set_visibility(view in ('temperature', 'pressure'))
        self.top.set_visibility(view in ('temperature', 'pressure'))
        self.log.set_visibility(view in ('temperature', 'pressure'))
        self.quantity.set_visibility(view != 'carbon')
        self.pressures.set_visibility(view in ('temperature', 'carbon'))
        self.temperature.set_visibility(view == 'pressure')
        self.item.set_visibility(view == 'map')

    def draw(self, _=None):
        self.controls()
        tables = self.tables
        if tables is None or tables.summary.empty:
            self.message.set_text('Choose a run.' if tables is None else 'No points solved yet.')
            self.plot.update_figure(plots.go.Figure())
            return
        solved = int((tables.summary.status == 'solved').sum()) if 'status' in tables.summary else len(tables.summary)
        self.message.set_text(f'{tables.label}: {solved} of {len(tables.summary)} points solved')
        view, quantity = self.view.value, self.quantity.value
        try:
            if view == 'temperature':
                comparison, label = self._comparison(tables, quantity)
                fig = plots.species_vs_temperature(tables.species, tables.summary, self.group.value, quantity,
                                                   self.pressures.value or None, int(self.top.value or 12),
                                                   self.log.value, comparison=comparison, comparison_label=label)
            elif view == 'pressure':
                fig = plots.species_vs_pressure(tables.species, tables.summary, self.group.value,
                                                self.temperature.value, quantity, int(self.top.value or 12),
                                                self.log.value)
            elif view == 'map':
                fig = plots.tp_map(tables, self.item.value or 'gas_carbon_percent', quantity)
            else:
                fig = plots.carbon_vs_temperature(tables.phases, tables.summary, self.pressures.value or None)
        except Exception as exc:
            self.message.set_text(f'Cannot draw this view: {exc}')
            fig = plots.go.Figure()
        self.plot.update_figure(fig)

    @staticmethod
    def _comparison(tables, quantity):
        settings = tables.case.comparison if tables.case is not None else None
        if tables.comparison is None or not settings:
            return None, ''
        if settings.get('quantity', 'mol_per_kg_dry_feed') != quantity:
            return None, ''
        return tables.comparison, settings.get('label') or 'Reference'

    def fill_tables(self, _=None):
        tables = self.tables
        empty = pd.DataFrame()
        show_table(self.summary_table, tables.summary if tables else empty)
        species = tables.species if tables else empty
        text = (self.species_filter.value or '').strip().casefold()
        if text and not species.empty:
            species = species[species.species.str.casefold().str.contains(text, regex=False)]
        show_table(self.species_table, species)
        show_table(self.phases_table, tables.phases if tables else empty)

    def download(self, name):
        if self.tables is None:
            ui.notify('Choose a run first.')
            return
        df = getattr(self.tables, name)
        stem = self.tables.case.name if self.tables.case is not None else 'results'
        ui.download.content(df.to_csv(index=False).encode(), f'{stem}_{name}.csv')
