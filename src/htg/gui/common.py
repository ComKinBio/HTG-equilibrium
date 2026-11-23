"""Small helpers shared by the tabs: tables and status text."""
import math
import time

from nicegui import ui


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
