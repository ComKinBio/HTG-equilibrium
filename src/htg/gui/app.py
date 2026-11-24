"""The app: ``htg gui`` (own window) or ``htg gui --browser``.

One shared state (model catalogue, worker pool, preferences, the case being
edited) and a page with five tabs: case editor, queue, results, compare, species.
"""
import os
from pathlib import Path

from nicegui import app, ui

from ..model import Model
from ..runner import WorkerPool
from ..settings import Settings
from .common import table_rows, worker_summary  # noqa: F401  (table_rows is used by tests)
from .compare_tab import CompareTab
from .document import CaseDocument
from .editor import CaseEditor
from .queue_tab import QueueTab
from .results_tab import ResultsTab
from .runs import live_tables, load_saved, output_folder, saved_runs
from .species_tab import SpeciesTab, add_species  # noqa: F401  (add_species is used by tests)
from .window import native_or_browser, quiet_qt, use_qt


class State:
    """Everything shared by the tabs (one local user).

    settings  app preferences (htg.settings.Settings): read once here, written by
              the tabs with settings.set(...). Default: Settings.load().
    """

    def __init__(self, model, pool, results_dir='results', cases_dir='cases', settings=None):
        self.pool = pool
        self.settings = settings if settings is not None else Settings.load()
        self.results_dir = Path(results_dir).resolve()
        self.cases_dir = Path(cases_dir).resolve()
        self.document = CaseDocument.new()
        self.use_model(model)

    def use_model(self, model):
        """Use ``model`` (for example after the user species changed): the species lists follow it."""
        self.model = model
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
    # The tabs are in the (fixed) header so they stay visible when a long tab is scrolled.
    with ui.header().classes('column items-stretch gap-0 p-0 bg-primary'):
        with ui.row().classes('items-center w-full px-4 pt-2'):
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


def run_app(results_dir='results', cases_dir='cases', workers=None, port=None, show=True, model_args=None,
            worker_target=None, warm_up=True, host='127.0.0.1', browser=False, settings=None):
    """Start the worker pool and serve the app until the window/terminal is closed.

    browser  False (default): open the app in its own window (pywebview); if that
             cannot start, print why and use the browser. True: browser tab.
    show     browser mode only: open a browser tab (False: only serve).
    port     default: 8080 in the browser, a free port for the window.
    host     '127.0.0.1' (default): only this computer can open the app. The app reads
             and writes files, so share it on a network only on purpose.
    """
    native = native_or_browser(browser or not show)
    model = Model.load(**(model_args or {}))
    pool = WorkerPool(workers, model.load_args, warm_up=warm_up, worker_target=worker_target)
    state = State(model, pool, results_dir, cases_dir, settings=settings)
    stopped = []

    def stop():
        if not stopped:
            stopped.append(True)
            pool.shutdown(wait=False)
            print('Stopped. Finished runs are saved in', state.results_dir)

    app.on_startup(pool.start)
    # Closing the window stops the server and runs on_shutdown, then NiceGUI exits the
    # process: code after ui.run() does not run in that case.
    app.on_shutdown(stop)
    print(f'Starting {pool.target} worker(s); each needs about a minute before the first point.')
    if native:
        if use_qt():
            app.native.start_args['gui'] = 'qt'
            quiet_qt()
        window = {'native': True, 'window_size': (1400, 900), 'port': port}
    else:
        window = {'show': show, 'port': port or 8080}
    try:
        ui.run(lambda: build_page(state), host=host, reload=False, title='HTG equilibrium',
               show_welcome_message=not native, **window)
    except KeyboardInterrupt:
        pass
    finally:
        stop()
