"""Compare tab: one quantity from several runs."""
import pandas as pd
from nicegui import ui

from ..constants import PHASES
from . import plots
from .common import data_table, show_table
from .runs import case_differences


class CompareTab:
    def __init__(self, state):
        self.state = state
        with ui.row().classes('w-full items-end gap-3'):
            self.runs = ui.select(state.sources(), label='Runs to compare', multiple=True, with_input=True,
                                  on_change=self.load).classes('w-[40rem]').props('use-chips')
            ui.button(icon='refresh', on_click=self.refresh_sources).props('flat round')
        with ui.row().classes('w-full items-end gap-3'):
            self.key = ui.select({}, label='Species or carbon share', with_input=True,
                                 on_change=self.draw).classes('w-80')
            self.quantity = ui.select(plots.QUANTITIES, label='Quantity', value='mol_per_kg_dry_feed',
                                      on_change=self.draw).classes('w-56')
            self.log = ui.switch('Log scale', value=False, on_change=self.draw)
            self.reference = ui.select({}, label='Comparison data', on_change=self.draw).classes('w-72')
        self.message = ui.label().classes('text-grey-7')
        self.plot = ui.plotly(plots.go.Figure()).classes('w-full h-[32rem]')
        ui.label('Settings that differ').classes('text-h6')
        self.differences = data_table(pd.DataFrame(), rows_per_page=20)
        self.loaded = []
        ui.timer(5.0, self.refresh_sources)

    def refresh_sources(self):
        sources = self.state.sources()
        if sources != self.runs.options:
            self.runs.set_options(sources, value=[v for v in self.runs.value or [] if v in sources])

    def load(self, _=None):
        self.loaded = []
        for source in self.runs.value or []:
            try:
                tables = self.state.tables(source)
            except Exception as exc:
                ui.notify(f'Could not read {source}: {exc}', type='negative')
                continue
            if tables is not None:
                self.loaded.append(tables)
        keys = {f'carbon:{p}': f'Carbon in {plots.PHASE_LABELS[p]}' for p in PHASES}
        for tables in self.loaded:
            keys.update({k: k.replace(':', ': ') for k in tables.species_keys})
        self.key.set_options(keys, value=self.key.value if self.key.value in keys else 'carbon:gas')
        references = {'': '(none)'}
        for tables in self.loaded:
            settings = tables.case.comparison if tables.case is not None else None
            if tables.comparison is not None and settings:
                references[tables.label] = f'{settings.get("label") or "Reference"} (from {tables.label})'
        self.reference.set_options(references, value=self.reference.value
                                   if self.reference.value in references else '')
        show_table(self.differences, case_differences({t.label: t.case for t in self.loaded}))
        self.draw()

    def draw(self, _=None):
        if not self.loaded:
            self.message.set_text('Choose two or more runs.')
            self.plot.update_figure(plots.go.Figure())
            return
        self.message.set_text(f'{len(self.loaded)} run(s)')
        comparison, label = None, 'Reference'
        source = next((t for t in self.loaded if t.label == self.reference.value), None)
        if source is not None:
            settings = source.case.comparison
            if settings.get('quantity', 'mol_per_kg_dry_feed') == self.quantity.value or \
                    (self.key.value or '').startswith('carbon:'):
                comparison, label = source.comparison, settings.get('label') or 'Reference'
            else:
                self.message.set_text('The comparison data are in another unit than the chosen quantity.')
        fig = plots.compare_runs(self.loaded, self.key.value or 'carbon:gas', self.quantity.value,
                                 self.log.value, comparison, label)
        self.plot.update_figure(fig)
