"""The browser app: ``htg gui``.

One shared state (model catalogue, worker pool, the case being edited) and a
page with five tabs: case editor, queue, results, compare, species.
"""
import math
import os
import time
from pathlib import Path

import pandas as pd
from nicegui import app, ui

from ..constants import PHASES
from ..model import Model
from ..runner import WorkerPool
from . import plots
from .document import CaseDocument
from .editor import PHASE_NAMES, CaseEditor
from .runs import case_differences, live_tables, load_saved, output_folder, saved_runs


class State:
    """Everything shared by the tabs (one local user)."""

    def __init__(self, model, pool, results_dir='results', cases_dir='cases'):
        self.model = model
        self.pool = pool
        self.results_dir = Path(results_dir).resolve()
        self.cases_dir = Path(cases_dir).resolve()
        self.document = CaseDocument.new()
        catalogue = model.species()
        molecular = catalogue[catalogue.Phase == 'molecular'].Species.tolist()
        ions = catalogue[catalogue.Phase == 'aqueous ion'].Species.tolist()
        solids = catalogue[catalogue.Phase == 'solid'].Species.tolist()
        self.phase_options = {'gas': molecular, 'aqueous': molecular + ions,
                              'organic_liquid': molecular, 'solid': solids}
        registry = model.registry
        self.feed_species = sorted({*registry.all_molecular, *registry.solid})

    def submit(self, case):
        return self.pool.submit(case, output_folder(self.results_dir, case.name))

    def comparison_files(self, base_dir):
        """CSV files below the cases folder: {path relative to base_dir: shown name}."""
        if not self.cases_dir.exists():
            return {}
        files = sorted(self.cases_dir.glob('**/*.csv'))
        return {os.path.relpath(p, base_dir).replace(os.sep, '/'): str(p.relative_to(self.cases_dir))
                for p in files}

    def sources(self):
        """Runs to show: this session's runs (newest first), then saved results folders."""
        live = {f'live:{run.id}': f'#{run.id} {run.name} ({run.status}, {run.done_points}/{run.total})'
                for run in reversed(self.pool.runs)}
        session_dirs = {run.output_dir.resolve() for run in self.pool.runs if run.output_dir}
        saved = {f'saved:{folder}': str(folder.relative_to(self.results_dir))
                 for folder in saved_runs(self.results_dir) if folder.resolve() not in session_dirs}
        return {**live, **saved}

    def run(self, source):
        if source and source.startswith('live:'):
            run_id = int(source[5:])
            return next((r for r in self.pool.runs if r.id == run_id), None)
        return None

    def tables(self, source):
        run = self.run(source)
        if run is not None:
            return live_tables(run)
        if source and source.startswith('saved:'):
            return load_saved(source[6:], self.results_dir)
        return None


# ---------------------------------------------------------------- helpers
def table_rows(df, digits=6):
    """DataFrame -> rows for ui.table: NaN as empty, floats rounded."""
    def tidy(value):
        if isinstance(value, float):
            return None if math.isnan(value) else float(f'{value:.{digits}g}')
        if hasattr(value, 'item'):  # numpy scalars
            return tidy(value.item())
        return value
    return [{'_row': i, **{k: tidy(v) for k, v in row.items()}} for i, row in enumerate(df.to_dict('records'))]


def table_columns(df):
    return [{'name': c, 'label': c, 'field': c, 'sortable': True, 'align': 'left'} for c in df.columns]


def data_table(df, rows_per_page=10, **kwargs):
    return ui.table(rows=table_rows(df), columns=table_columns(df), row_key='_row',
                    pagination=rows_per_page, **kwargs).classes('w-full').props('dense flat bordered')


def show_table(table, df):
    table.columns = table_columns(df)
    table.rows = table_rows(df)
    table.update()


def elapsed(run):
    if run.started is None:
        return ''
    seconds = (run.finished or time.time()) - run.started
    return f'{seconds / 60:.1f} min' if seconds >= 120 else f'{seconds:.0f} s'


def worker_summary(status):
    counts = {}
    for worker in status['workers']:
        counts[worker['state']] = counts.get(worker['state'], 0) + 1
    order = ['busy', 'idle', 'warming', 'starting', 'stopping', 'crashed']
    parts = [f'{counts[s]} {s}' for s in order if counts.get(s)]
    text = 'Workers: ' + (', '.join(parts) or 'none')
    if status['queued_pieces']:
        text += f' · {status["queued_pieces"]} piece(s) waiting'
    return text


# ---------------------------------------------------------------- queue tab
class RunCard:
    """Progress of one run in the queue tab."""

    def __init__(self, run, pool, show_results):
        self.run = run
        with ui.card().classes('w-full') as self.card:
            with ui.row().classes('w-full items-center'):
                ui.label(f'#{run.id} {run.name}').classes('text-subtitle1 text-weight-medium')
                self.status = ui.badge(run.status)
                ui.space()
                self.cancel = ui.button('Cancel', icon='stop', on_click=lambda: pool.cancel(run)) \
                    .props('flat dense color=negative')
                ui.button('Results', icon='insights', on_click=lambda: show_results(f'live:{run.id}')) \
                    .props('flat dense')
            self.progress = ui.linear_progress(value=0, show_value=False).props('size=10px')
            self.text = ui.label().classes('text-caption')
            with ui.expansion('Log').classes('w-full'):
                self.log = ui.log(max_lines=500).classes('w-full h-40')
        self.logged = 0
        self.update()

    def update(self):
        run = self.run
        self.status.set_text(run.status)
        self.status.props(f'color={dict(done="positive", cancelled="grey", running="primary").get(run.status, "orange")}')
        self.progress.set_value(run.done_points / run.total if run.total else 1)
        text = f'{run.done_points} of {run.total} points'
        if run.failed_points:
            text += f', {run.failed_points} failed'
        if run.started:
            text += f' · {elapsed(run)}'
        if run.finished_running and run.output_dir:
            text += f' · saved to {run.output_dir}'
        self.text.set_text(text)
        self.cancel.set_visibility(not run.finished_running)
        for line in run.log[self.logged:]:
            self.log.push(line)
        self.logged = len(run.log)


class QueueTab:
    def __init__(self, state, show_results):
        self.state = state
        self.show_results = show_results
        self.cards = {}
        pool = state.pool
        with ui.row().classes('w-full items-end gap-4'):
            self.workers = ui.number('Parallel workers', value=pool.target, min=1, max=os.cpu_count() or 1,
                                     step=1, precision=0).classes('w-40')
            ui.button('Apply', on_click=self._apply).props('outline')
            ui.label('Each worker uses about 1.5 GB of memory and needs about a minute to start. '
                     'Independent pieces (single points, or one continuation chain per pressure) '
                     'are handed to free workers in the order they were queued.') \
                .classes('text-caption text-grey-7 grow')
        self.worker_table = ui.table(rows=[], columns=[
            {'name': k, 'label': label, 'field': k, 'align': 'left'}
            for k, label in (('id', 'Worker'), ('state', 'State'), ('case', 'Case'), ('message', 'Message'))],
            row_key='id').classes('w-full').props('dense flat bordered')
        ui.label('Runs').classes('text-h6 q-mt-md')
        self.empty = ui.label('Nothing queued yet: open or edit a case and press "Add to queue".') \
            .classes('text-grey-7')
        self.runs_box = ui.column().classes('w-full gap-2')
        self.update()
        ui.timer(1.0, self.update)

    def _apply(self):
        if self.workers.value:
            self.state.pool.set_workers(int(self.workers.value))
            ui.notify(f'Using {int(self.workers.value)} worker(s)')

    def update(self):
        status = self.state.pool.status()
        self.worker_table.rows = sorted(status['workers'], key=lambda w: w['id'])
        self.worker_table.update()
        for run in self.state.pool.runs:
            if run.id not in self.cards:
                with self.runs_box:
                    card = RunCard(run, self.state.pool, self.show_results)
                card.card.move(self.runs_box, target_index=0)  # newest first
                self.cards[run.id] = card
            elif not (self.cards[run.id].run.finished_running and self.cards[run.id].status.text == run.status):
                self.cards[run.id].update()
        self.empty.set_visibility(not self.cards)


# ---------------------------------------------------------------- results tab
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


# ---------------------------------------------------------------- compare tab
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


# ---------------------------------------------------------------- species tab
TARGETS = {'molecular': ['gas', 'aqueous', 'organic_liquid'], 'aqueous ion': ['aqueous'], 'solid': ['solid']}


class SpeciesTab:
    def __init__(self, state, editor):
        self.state = state
        self.editor = editor
        catalogue = state.model.species(include_unavailable=True)
        self.catalogue = catalogue[['Species', 'Formula', 'Phase', 'Ready', 'Reason']].copy()
        ui.label('Species with complete data are "ready". Select one to add it to the case being edited.') \
            .classes('text-caption text-grey-7')
        with ui.row().classes('w-full items-end gap-3'):
            self.search = ui.input('Search name or formula', on_change=self.fill).classes('w-64')
            self.phase = ui.select(['all', *sorted(self.catalogue.Phase.unique())], value='all', label='Kind',
                                   on_change=self.fill).classes('w-56')
            self.unavailable = ui.switch('Show unavailable', value=False, on_change=self.fill)
        with ui.row().classes('items-center gap-2'):
            self.selected_label = ui.label('Nothing selected').classes('text-weight-medium')
            self.buttons = {}
            for phase, label in PHASE_NAMES.items():
                self.buttons[phase] = ui.button(f'Add to {label.lower()}', icon='add',
                                                on_click=lambda p=phase: self.add(p)).props('outline dense')
            self.buttons['feed'] = ui.button('Add to feed species', icon='add',
                                             on_click=lambda: self.add('feed')).props('outline dense')
        self.table = ui.table(rows=[], columns=table_columns(self.catalogue), row_key='key', selection='single',
                              pagination=20, on_select=self.selected).classes('w-full').props('dense flat bordered')
        self.choice = None
        self.fill()
        self.selected()

    def fill(self, _=None):
        df = self.catalogue
        if not self.unavailable.value:
            df = df[df.Ready]
        if self.phase.value and self.phase.value != 'all':
            df = df[df.Phase == self.phase.value]
        text = (self.search.value or '').strip().casefold()
        if text:
            df = df[df.Species.astype(str).str.casefold().str.contains(text, regex=False)
                    | df.Formula.astype(str).str.casefold().str.contains(text, regex=False)]
        rows = table_rows(df)
        for row in rows:
            row['key'] = f'{row["Phase"]}:{row["Species"]}'
        self.table.rows = rows
        self.table.update()

    def selected(self, _=None):
        rows = self.table.selected
        self.choice = rows[0] if rows else None
        allowed = set(TARGETS.get(self.choice['Phase'], [])) if self.choice and self.choice['Ready'] else set()
        for phase, button in self.buttons.items():
            enabled = (self.choice is not None and self.choice['Ready']
                       and self.choice['Species'] in self.state.feed_species) if phase == 'feed' \
                else phase in allowed
            button.set_enabled(enabled)
        if self.choice is None:
            self.selected_label.set_text('Nothing selected')
        else:
            note = '' if self.choice['Ready'] else f' — not usable: {self.choice["Reason"]}'
            self.selected_label.set_text(f'{self.choice["Species"]} ({self.choice["Formula"]}){note}')

    def add(self, target):
        if self.choice is None:
            return
        name = self.choice['Species']
        try:
            where = add_species(self.state.document.data, target, name)
        except ValueError as exc:
            ui.notify(str(exc))
            return
        self.editor.form.refresh()
        self.editor.changed()
        ui.notify(f'Added {name} to {where}', type='positive')


def add_species(data, target, name):
    """Add ``name`` to a phase list ('gas', ...) or to the feed species ('feed') of case data."""
    if target == 'feed':
        feed = data.setdefault('feed', {})
        feed.pop('elements_wt_percent', None)
        feed.setdefault('species_moles', {}).setdefault(name, 1.0)
        return 'the feed'
    phases = data.setdefault('phases', {})
    current = list(phases.get(target) or [])
    if name in current:
        raise ValueError(f'{name} is already in the {PHASE_NAMES[target].lower()} phase.')
    phases[target] = current + [name]
    return f'the {PHASE_NAMES[target].lower()} phase'


# ---------------------------------------------------------------- page
def build_page(state):
    ui.page_title('HTG equilibrium')
    with ui.header().classes('items-center bg-primary'):
        ui.label('HTG equilibrium').classes('text-h6')
        ui.space()
        status = ui.label().classes('text-body2')
        ui.timer(1.0, lambda: status.set_text(worker_summary(state.pool.status())))
    with ui.tabs().classes('w-full') as tabs:
        editor_tab = ui.tab('Case', icon='edit_note')
        queue_tab = ui.tab('Queue', icon='queue')
        results_tab = ui.tab('Results', icon='insights')
        compare_tab = ui.tab('Compare', icon='compare_arrows')
        species_tab = ui.tab('Species', icon='science')
    results = None

    def show_results(source):
        tabs.set_value(results_tab)
        results.select(source)

    with ui.tab_panels(tabs, value=editor_tab).classes('w-full max-w-screen-2xl mx-auto'):
        with ui.tab_panel(editor_tab):
            editor = CaseEditor(state, on_queued=lambda run: tabs.set_value(queue_tab))
        with ui.tab_panel(queue_tab):
            QueueTab(state, show_results)
        with ui.tab_panel(results_tab):
            results = ResultsTab(state)
        with ui.tab_panel(compare_tab):
            CompareTab(state)
        with ui.tab_panel(species_tab):
            SpeciesTab(state, editor)


def run_app(results_dir='results', cases_dir='cases', workers=None, port=8080, show=True, model_args=None,
            worker_target=None, warm_up=True, host='127.0.0.1'):
    """Start the worker pool and serve the app until the window/terminal is closed.

    host  '127.0.0.1' (default): only this computer can open the app. The app reads
          and writes files, so share it on a network only on purpose.
    """
    model = Model.load(**(model_args or {}))
    pool = WorkerPool(workers, model.load_args, warm_up=warm_up, worker_target=worker_target)
    state = State(model, pool, results_dir, cases_dir)
    app.on_startup(pool.start)
    app.on_shutdown(lambda: pool.shutdown(wait=False))
    print(f'Starting {pool.target} worker(s); each needs about a minute before the first point.')
    try:
        ui.run(lambda: build_page(state), host=host, port=port, show=show, reload=False, title='HTG equilibrium',
               show_welcome_message=True)
    except KeyboardInterrupt:
        pass
    finally:
        pool.shutdown(wait=False)
        print('Stopped. Finished runs are saved in', state.results_dir)
