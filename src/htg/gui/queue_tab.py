"""Queue tab: worker pool and run progress."""
import os

from nicegui import ui

from .common import elapsed


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
