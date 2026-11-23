"""Species tab: browse the species catalogue and add species to the case."""
from nicegui import ui

from .common import table_columns, table_rows
from .editor import PHASE_NAMES


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
