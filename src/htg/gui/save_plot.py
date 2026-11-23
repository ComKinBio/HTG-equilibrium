"""The "Save plot…" dialog: layout, format, size and a preview.

The file is drawn with Matplotlib (htg.plotrender.render_matplotlib) from the same
PlotSpec as the plot on screen. Width, height, resolution, format and the legend
choice are remembered per layout in the app preferences.
"""
import base64
import io

from nicegui import ui

from .. import plotspec as ps
from ..plotrender import render_matplotlib
from .common import save_bytes

PREVIEW_PX = 440  # width of the preview image
SIZE_LIMITS_MM = (10, 1000)
DPI_LIMITS = (50, 2400)


def render_bytes(spec, layout, width_mm, height_mm, fmt, dpi=ps.DEFAULT_DPI, shared_legend=True,
                 theme=ps.DEFAULT_THEME, line_width=ps.DEFAULT_LINE_WIDTH, panel=0):
    """The plot file as bytes (fmt 'PNG', 'SVG' or 'PDF')."""
    buffer = io.BytesIO()
    render_matplotlib(spec, buffer, layout, width_mm, height_mm, fmt, dpi=dpi, shared_legend=shared_legend,
                      theme=theme, line_width=line_width, panel=panel)
    return buffer.getvalue()


class SavePlotDialog:
    """Dialog for one PlotSpec.

    settings   htg.settings.Settings (export defaults per layout)
    spec       the PlotSpec to save
    theme, line_width   the look chosen on the tab
    stem       file name without suffix
    layout     layout to start with (default: as shown on the tab)
    """

    def __init__(self, settings, spec, theme, line_width, stem='plot', layout=None):
        self.settings = settings
        self.spec = spec
        self.theme = theme
        self.line_width = line_width
        self.stem = stem
        self.layouts = {option.key: option for option in ps.layout_options(spec)}
        if layout not in self.layouts:
            layout = 'grid' if 'grid' in self.layouts and spec.rows > 1 else 'one'
        self._quiet = True
        panels = {i: spec.panels[i].title or spec.title for i in self.layouts['one'].panels
                  if not spec.panels[i].empty}
        with ui.dialog() as self.dialog, ui.card().classes('w-[56rem] max-w-full'):
            ui.label('Save plot').classes('text-h6')
            with ui.row().classes('w-full no-wrap gap-6'):
                with ui.column().classes('w-80 shrink-0 gap-2'):
                    ui.label('Layout').classes('text-weight-medium')
                    self.layout = ui.radio({k: o.label for k, o in self.layouts.items()}, value=layout,
                                           on_change=self._layout_changed).props('dense')
                    self.caption = ui.label().classes('text-caption text-grey-7')
                    self.panel = ui.select(panels, label='Panel', value=next(iter(panels), 0),
                                           on_change=self._changed).classes('w-full')
                    ui.label('Format').classes('text-weight-medium')
                    self.format = ui.toggle(list(ps.FORMATS), value='PNG', on_change=self._changed)
                    with ui.row().classes('w-full no-wrap gap-2'):
                        self.width = ui.number('Width [mm]', min=SIZE_LIMITS_MM[0], max=SIZE_LIMITS_MM[1],
                                               on_change=self._changed).classes('w-32')
                        self.height = ui.number('Height [mm]', min=SIZE_LIMITS_MM[0], max=SIZE_LIMITS_MM[1],
                                                on_change=self._changed).classes('w-32')
                    self.dpi = ui.number('Resolution [dpi]', min=DPI_LIMITS[0], max=DPI_LIMITS[1], step=50,
                                         precision=0, on_change=self._changed).classes('w-32')
                    self.shared = ui.checkbox('One shared legend', value=True, on_change=self._changed)
                with ui.column().classes('grow items-center gap-2'):
                    self.preview = ui.image().classes('w-full border').style(f'max-width:{PREVIEW_PX}px')
                    self.summary = ui.label().classes('text-body2')
                    self.problem = ui.label().classes('text-negative text-caption')
            with ui.row().classes('w-full justify-end gap-2'):
                ui.button('Cancel', on_click=self.dialog.close).props('flat no-caps')
                self.save_button = ui.button('Save…', icon='save', on_click=self.save).props('no-caps')
        self._load_defaults(layout)
        self._quiet = False
        self._changed()

    # ------------------------------------------------------------ values
    def _load_defaults(self, layout):
        values = self.settings.export_defaults(layout)
        self.width.set_value(values['width_mm'])
        self.height.set_value(values['height_mm'])
        self.dpi.set_value(values['dpi'])
        self.format.set_value(values['format'].upper())
        self.shared.set_value(values['shared_legend'])

    def values(self):
        """The current choices, or a text saying which field is wrong."""
        layout = self.layout.value
        width, height, dpi = self.width.value, self.height.value, self.dpi.value
        low, high = SIZE_LIMITS_MM
        for name, value in (('Width', width), ('Height', height)):
            if value is None or not low <= value <= high:
                return f'{name}: enter a value from {low} to {high} mm.'
        if self.format.value == 'PNG' and (dpi is None or not DPI_LIMITS[0] <= dpi <= DPI_LIMITS[1]):
            return f'Resolution: enter a value from {DPI_LIMITS[0]} to {DPI_LIMITS[1]} dpi.'
        return {'layout': layout, 'width_mm': float(width), 'height_mm': float(height),
                'dpi': int(dpi or ps.DEFAULT_DPI), 'fmt': self.format.value, 'shared_legend': bool(self.shared.value),
                'panel': int(self.panel.value or 0)}

    def file_bytes(self, values=None, dpi=None):
        values = values or self.values()
        return render_bytes(self.spec, values['layout'], values['width_mm'], values['height_mm'],
                            values['fmt'] if dpi is None else 'PNG', dpi or values['dpi'], values['shared_legend'],
                            self.theme, self.line_width, values['panel'])

    # ------------------------------------------------------------ events
    def _layout_changed(self, _=None):
        if self._quiet:
            return
        self._quiet = True
        self._load_defaults(self.layout.value)
        self._quiet = False
        self._changed()

    def _changed(self, _=None):
        if self._quiet:
            return
        option = self.layouts[self.layout.value]
        self.caption.set_text(option.caption)
        self.caption.set_visibility(bool(option.caption))
        self.panel.set_visibility(self.layout.value == 'one' and len(self.panel.options) > 1)
        self.dpi.set_visibility(self.format.value == 'PNG')
        values = self.values()
        if isinstance(values, str):
            self.problem.set_text(values)
            self.save_button.disable()
            return
        self.problem.set_text('')
        self.save_button.enable()
        self.summary.set_text(ps.summary_line(self.spec, values['layout'], values['width_mm'], values['height_mm'],
                                              values['fmt'], values['panel']))
        dpi = max(20, min(200, PREVIEW_PX * 25.4 / values['width_mm']))
        try:
            data = self.file_bytes(values, dpi=dpi)
        except Exception as exc:  # the preview must never break the dialog
            self.problem.set_text(f'No preview: {exc}')
            return
        self.preview.set_source('data:image/png;base64,' + base64.b64encode(data).decode())

    async def save(self):
        values = self.values()
        if isinstance(values, str):
            ui.notify(values, type='warning')
            return
        self.settings.set_export_defaults(values['layout'], width_mm=values['width_mm'],
                                          height_mm=values['height_mm'], dpi=values['dpi'],
                                          format=values['fmt'].lower(), shared_legend=values['shared_legend'])
        try:
            data = self.file_bytes(values)
        except Exception as exc:
            ui.notify(f'Could not draw the plot: {exc}', type='negative', multi_line=True)
            return
        where = await save_bytes(data, f'{self.stem}.{values["fmt"].lower()}')
        if where is None:
            return
        self.dialog.close()
        ui.notify(f'Saved {where}')

    def open(self):
        self.dialog.open()
        return self
