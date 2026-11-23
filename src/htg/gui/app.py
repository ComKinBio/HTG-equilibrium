"""The browser app: ``htg gui``.

One shared state (model catalogue, worker pool, the case being edited) and a
page with five tabs: case editor, queue, results, compare, species.
"""
import os
from pathlib import Path

from nicegui import app, ui

from ..model import Model
from ..runner import WorkerPool
from .common import table_rows, worker_summary  # noqa: F401  (table_rows is used by tests)
from .compare_tab import CompareTab
from .document import CaseDocument
from .editor import CaseEditor
from .queue_tab import QueueTab
from .results_tab import ResultsTab
from .runs import live_tables, load_saved, output_folder, saved_runs
from .species_tab import SpeciesTab, add_species  # noqa: F401  (add_species is used by tests)


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
