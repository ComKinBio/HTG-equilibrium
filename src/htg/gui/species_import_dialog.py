"""The "Import from file…" dialog of the Species tab (see htg.species_import).

    ImportDialog(model, on_imported).open()

The file comes from a native file dialog in the app window, else from an upload
field in the browser. on_imported(count) is called after the user database
changed; it may be a coroutine function.
"""
import inspect
from pathlib import Path

from nicegui import app, ui
from nicegui import run as background

from ..species_import import FORMATS, NOTE, apply_import, detect_format, preview_import, read_text

COLUMNS = [{'name': c, 'label': c, 'field': c, 'align': 'left'}
           for c in ('Species', 'Formula', 'What happens', 'After import')]
SELECTION_SLOT = ('<q-td auto-width><q-checkbox v-model="props.selected" dense '
                  ':disable="!props.row.selectable" /></q-td>')


class ImportDialog:
    def __init__(self, model, on_imported):
        self.model = model
        self.on_imported = on_imported
        self.text = None
        self.rows = []
        self._loading = False
        self.dialog = ui.dialog()
        with self.dialog, ui.card().classes('w-[68rem] max-w-full'):
            ui.label('Import from file').classes('text-h6')
            with ui.row().classes('w-full items-center gap-4'):
                if app.native.main_window is not None:
                    ui.button('Choose file…', on_click=self.pick_native).props('outline no-caps')
                else:
                    ui.upload(label='Choose file', auto_upload=True, max_files=1, on_upload=self.uploaded) \
                        .props('flat bordered').classes('w-72')
                with ui.column().classes('gap-0'):
                    self.file_label = ui.label('No file chosen').classes('text-body2')
                    self.format_note = ui.label().classes('text-caption text-grey-7')
                ui.space()
                self.format = ui.select(dict(FORMATS), label='Format', on_change=self.format_changed) \
                    .props('dense').classes('w-56')
            self.message = ui.label().classes('text-negative')
            self.table = ui.table(rows=[], columns=COLUMNS, row_key='key', selection='multiple',
                                  pagination=15, on_select=self.selection_changed) \
                .classes('w-full').props('dense flat bordered no-data-label="No species to show"')
            self.table.add_slot('body-selection', SELECTION_SLOT)
            ui.label(NOTE).classes('text-caption text-grey-7')
            with ui.row().classes('w-full'):
                ui.space()
                ui.button('Cancel', on_click=self.dialog.close).props('flat no-caps')
                self.import_button = ui.button('Import 0 species', on_click=self.do_import).props('no-caps')
        self.update_button()

    def open(self):
        self.dialog.open()
        return self

    # ------------------------------------------------------------ the file
    async def pick_native(self):
        paths = await app.native.main_window.create_file_dialog(allow_multiple=False)
        if not paths:
            return
        path = paths[0] if isinstance(paths, (list, tuple)) else paths
        try:
            text = read_text(path)
        except OSError as exc:
            self.message.set_text(f'Could not read {path}: {exc}')
            return
        await self.load_text(text, Path(path).name)

    async def uploaded(self, e):
        data = await e.file.read()
        await self.load_text(data.decode('utf-8-sig', errors='replace'), e.file.name)

    async def load_text(self, text, name):
        """Show the preview for file content ``text`` (format detected from the content)."""
        self.text = text
        self.file_label.set_text(name)
        fmt = detect_format(text)
        self.format_note.set_text('Format detected from the content' if fmt else
                                  'Format not recognised. Choose the format.')
        self._loading = True
        try:
            self.format.set_value(fmt)
        finally:
            self._loading = False
        await self.preview()

    async def format_changed(self, _):
        if not self._loading:
            await self.preview()

    async def preview(self):
        self.rows = []
        self.message.set_text('')
        fmt = self.format.value
        if self.text is not None and fmt:
            try:
                self.rows = await background.io_bound(preview_import, self.text, self.model, fmt)
            except Exception as exc:  # a file in another format can fail in many ways
                self.message.set_text(f'Could not read the file as {FORMATS[fmt]}: {exc}')
            if not self.rows and not self.message.text:
                self.message.set_text('No species found. Check the format.')
        table_rows = [row.as_row() for row in self.rows]
        self.table.rows = table_rows
        self.table.selected = [row for row in table_rows if row['selected']]
        self.table.update()
        self.update_button()

    # ------------------------------------------------------------ ticks and import
    def selection_changed(self, _=None):
        allowed = [row for row in self.table.selected if row.get('selectable')]
        if len(allowed) != len(self.table.selected):
            self.table.selected = allowed
            self.table.update()
        self.update_button()

    def chosen(self):
        keys = {row['key'] for row in self.table.selected}
        return [row for row in self.rows if row.key in keys and row.selectable]

    def update_button(self):
        n = len(self.chosen())
        self.import_button.set_text(f'Import {n} species')
        self.import_button.set_enabled(n > 0)

    async def do_import(self):
        rows = self.chosen()
        if not rows:
            return
        try:
            count = apply_import(rows, self.model.user_data)
        except (OSError, ValueError) as exc:
            ui.notify(f'Not imported: {exc}', type='negative')
            return
        self.dialog.close()
        result = self.on_imported(count)
        if inspect.isawaitable(result):
            await result
