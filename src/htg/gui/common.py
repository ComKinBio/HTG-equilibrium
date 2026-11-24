"""Small helpers shared by the tabs: tables, status text, swatches and saving files."""
import math
import time
from pathlib import Path

from nicegui import app, ui


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


class OptionRadio(ui.radio):
    """A radio group whose options can be disabled one by one.

    set_disabled({value: reason}) disables those options; the reasons become the
    group's tooltip.
    """

    def __init__(self, options, *, value=None, on_change=None):
        self.disabled = {}
        super().__init__(options, value=value, on_change=on_change)
        with self:
            self.reason = ui.tooltip('')
        self.reason.set_visibility(False)

    def _update_options(self):
        super()._update_options()
        for option, value in zip(self._props['options'], self._values):
            if value in getattr(self, 'disabled', {}):
                option['disable'] = True

    def set_disabled(self, disabled):
        disabled = dict(disabled)
        if disabled == self.disabled:
            return
        self.disabled = disabled
        self._update_options()
        self.update()
        text = ' '.join(dict.fromkeys(disabled.values()))
        self.reason.set_text(text)
        self.reason.set_visibility(bool(text))


def swatch(colour, dash='solid', width=22):
    """A short line in a series colour (dashed or dotted when the colour repeats)."""
    style = {'solid': 'solid', 'dash': 'dashed', 'dot': 'dotted', 'dashdot': 'dashed'}.get(dash, 'solid')
    return ui.element('span').style(f'display:inline-block;width:{width}px;height:0;'
                                    f'border-top:3px {style} {colour};vertical-align:middle')


def palette_strip(container, colours):
    """Fill ``container`` with one small square per colour."""
    container.clear()
    with container:
        for colour in colours:
            ui.element('span').style(f'display:inline-block;width:12px;height:12px;background:{colour}')


FILE_TYPES = {'png': 'PNG image (*.png)', 'svg': 'SVG image (*.svg)', 'pdf': 'PDF document (*.pdf)',
              'zip': 'ZIP archive (*.zip)', 'csv': 'CSV table (*.csv)'}


async def save_bytes(data, filename):
    """Save ``data`` with a save dialog in the app window, or as a download in the browser.

    Returns the path written (window), the file name (download), or None if cancelled.
    """
    window = getattr(app.native, 'main_window', None)
    if window is None:
        ui.download.content(data, filename)
        return filename
    import webview

    kind = webview.FileDialog.SAVE if hasattr(webview, 'FileDialog') else webview.SAVE_DIALOG
    suffix = Path(filename).suffix.lstrip('.').lower()
    types = (FILE_TYPES[suffix],) if suffix in FILE_TYPES else ()
    chosen = await window.create_file_dialog(kind, save_filename=filename, file_types=types)
    if not chosen:
        return None
    path = chosen if isinstance(chosen, str) else chosen[0]
    Path(path).write_bytes(data)
    return path
