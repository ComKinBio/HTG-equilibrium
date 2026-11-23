"""Queue tab: worker pool and run progress."""
import os
from pathlib import Path

from nicegui import ui

from .. import units
from .common import elapsed

WORKERS_CAPTION = ('Number of points solved at the same time. Each worker needs about 1.5 GB of memory '
                   'and takes about a minute to start.')
EMPTY = 'No runs yet. Add a case to the queue from the Case tab.'
WORKER_COLUMNS = (('id', 'Worker'), ('state', 'State'), ('case', 'Case'), ('working_on', 'Working on'))
STATUS_COLOURS = {'done': 'positive', 'cancelled': 'grey', 'running': 'primary'}


def point_text(point, t_unit='C', p_unit='bar'):
    """'T = 400 °C, P = 240 bar' in the display units."""
    T, P = point
    return (f'T = {units.display(T, t_unit)} {units.label(t_unit)}, '
            f'P = {units.display(P, p_unit)} {units.label(p_unit)}')


def worker_rows(status, t_unit='C', p_unit='bar'):
    """Rows of the worker table: Worker · State · Case · Working on."""
    rows = []
    for worker in sorted(status['workers'], key=lambda w: w['id']):
        point = worker.get('point')
        working_on = point_text(point, t_unit, p_unit) if point is not None else worker.get('message', '')
        rows.append({'id': worker['id'], 'state': worker['state'], 'case': worker.get('case', ''),
                     'working_on': working_on})
    return rows


def saved_place(output_dir, results_dir=None):
    """'results/manure-20260102-030405' (relative to the folder above the results folder)."""
    path = Path(output_dir)
    if results_dir is not None:
        try:
            return str(path.relative_to(Path(results_dir).parent))
        except ValueError:
            pass
    return str(path)


def run_text(run, results_dir=None):
    """'8 of 13 points · 2.4 min' or '13 of 13 points, 1 failed · 4.0 min · saved in results/…'."""
    text = f'{run.done_points} of {run.total} points'
    if run.failed_points:
        text += f', {run.failed_points} failed'
    if run.started:
        text += f' · {elapsed(run)}'
    if run.finished_running and run.output_dir and (Path(run.output_dir) / 'summary.csv').exists():
        text += f' · saved in {saved_place(run.output_dir, results_dir)}'
    return text


class RunCard:
    """Progress of one run in the queue tab."""

    def __init__(self, run, pool, show_results, results_dir=None):
        self.run = run
        self.results_dir = results_dir
        with ui.card().classes('w-full') as self.card:
            with ui.row().classes('w-full items-center'):
                ui.label(f'#{run.id} {run.name}').classes('text-subtitle1 text-weight-medium')
                self.status = ui.badge(run.status)
                ui.space()
                self.cancel = ui.button('Cancel', icon='stop', on_click=lambda: pool.cancel(run)) \
                    .props('flat dense no-caps color=negative')
                ui.button('Results', icon='insights', on_click=lambda: show_results(f'live:{run.id}')) \
                    .props('flat dense no-caps')
            self.progress = ui.linear_progress(value=0, show_value=False).props('size=10px')
            self.text = ui.label().classes('text-caption')
            with ui.expansion('Log').classes('w-full'):
                self.log = ui.log(max_lines=500).classes('w-full h-40')
        self.logged = 0
        self.update()

    def update(self):
        run = self.run
        self.status.set_text(run.status)
        self.status.props(f'color={STATUS_COLOURS.get(run.status, "orange")}')
        self.progress.set_value(run.done_points / run.total if run.total else 1)
        self.text.set_text(run_text(run, self.results_dir))
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
            self.workers = ui.number('Workers', value=pool.target, min=1, max=os.cpu_count() or 1,
                                     step=1, precision=0).classes('w-32')
            ui.button('Apply', on_click=self._apply).props('outline no-caps')
            ui.label(WORKERS_CAPTION).classes('text-caption text-grey-7 grow')
        self.worker_table = ui.table(rows=[], columns=[
            {'name': k, 'label': label, 'field': k, 'align': 'left'} for k, label in WORKER_COLUMNS],
            row_key='id').classes('w-full').props('dense flat bordered')
        ui.label('Runs').classes('text-h6 q-mt-md')
        self.empty = ui.label(EMPTY).classes('text-grey-7')
        self.runs_box = ui.column().classes('w-full gap-2')
        self.update()
        ui.timer(1.0, self.update)

    def _apply(self):
        if self.workers.value:
            n = int(self.workers.value)
            self.state.pool.set_workers(n)
            ui.notify(f'Using {n} worker{"s" if n != 1 else ""}')

    def update(self):
        settings = self.state.settings
        rows = worker_rows(self.state.pool.status(), settings.get('temperature_unit'), settings.get('pressure_unit'))
        if rows != self.worker_table.rows:
            self.worker_table.rows = rows
            self.worker_table.update()
        results_dir = getattr(self.state, 'results_dir', None)
        for run in self.state.pool.runs:
            if run.id not in self.cards:
                with self.runs_box:
                    card = RunCard(run, self.state.pool, self.show_results, results_dir)
                card.card.move(self.runs_box, target_index=0)  # newest first
                self.cards[run.id] = card
            elif not (self.cards[run.id].run.finished_running and self.cards[run.id].status.text == run.status):
                self.cards[run.id].update()
        self.empty.set_visibility(not self.cards)
