"""Species tab: browse the species, add them to the case, and edit your own species database.

Saving, importing or removing species reloads the model with the same load
arguments, refreshes the species lists of the app and restarts the worker pool
so the workers use the new data (see WorkerPool.restart).
"""
import re

from nicegui import run as background
from nicegui import ui

from ..feed_v2 import elemental_masses
from ..formulas import parse_element_list
from ..model import Model
from ..species_catalogue import (PHASE_LABELS, SHOW_OPTIONS, SOURCE_OPTIONS, build_catalogue, count_text,
                                 filter_catalogue, table_rows)
from ..user_species import UserDatabase, spec_from_model
from .case_form import tidy
from .editor import PHASE_NAMES
from .species_dialog import SpeciesDialog
from .species_import_dialog import ImportDialog

ADD_NOTE = 'Adds it to the case open on the Case tab'
DATABASE_NOTE = ('New and edited species are saved in your own database; the shipped databases are '
                 'never changed.')
DATABASE_OFF = 'Your own species database is off (--no-user-species).'
RESTART_LATER = 'Runs already queued use the old species data. Later runs use the new data.'
COLUMNS = [
    {'name': 'species', 'label': 'Species', 'field': 'Species', 'sortable': True, 'align': 'left'},
    {'name': 'formula', 'label': 'Formula', 'field': 'Formula', 'sortable': True, 'align': 'left'},
    {'name': 'can', 'label': 'Can be in', 'field': 'Can be in', 'align': 'left'},
    {'name': 'data', 'label': 'Data', 'field': 'Data', 'align': 'left'},
    {'name': 'status', 'label': 'Status', 'field': 'Status', 'sortable': True, 'align': 'left'},
]
CAN_BE_IN_SLOT = ('<q-td :props="props"><q-chip v-for="p in props.value" :key="p" dense square size="sm" '
                  'color="blue-1" text-color="blue-10" :label="p" /></q-td>')
DATA_SLOT = ('<q-td :props="props"><template v-for="d in props.value" :key="d">'
             '<q-badge v-if="d === \'yours\'" color="primary" label="yours" class="q-mr-xs" />'
             '<span v-else class="q-mr-sm">{{ d }}</span></template></q-td>')
STATUS_SLOT = ('<q-td :props="props" :class="props.row.Usable ? \'\' : \'text-grey-7\'">'
               '{{ props.value }}</q-td>')
ROW_FIELDS = ('Species', 'Formula', 'Can be in', 'Data', 'Status', 'Usable', 'Feed', 'Yours', 'Overrides',
              'Kind')


def element_error(text):
    """Validation message for the "Contains elements" field (None when fine)."""
    try:
        parse_element_list(text)
    except ValueError as exc:
        return f'{exc} Use symbols such as C, H, O.'
    return None


def add_species(data, target, name):
    """Add ``name`` to a phase list ('gas', ...) or to the feed species ('feed') of case data.

    The feed part follows the case file format and changes with it (version 2).
    """
    if target == 'feed':
        return add_feed_species(data, name)
    phases = data.setdefault('phases', {})
    current = list(phases.get(target) or [])
    if name in current:
        raise ValueError(f'{name} is already in the {PHASE_NAMES[target].lower()} phase.')
    phases[target] = current + [name]
    return f'the {PHASE_NAMES[target].lower()} phase'


def add_feed_species(data, name):
    """Add ``name`` to the feed species of case data ([feed.species] of case-file format 2).

    An elemental-analysis feed becomes a species feed: its water stays (water is
    always present) and the first species gets the mass of the solids. Later
    species start at 1 kg. Salts are kept.
    """
    feed = data.setdefault('feed', {})
    species = feed.get('species')
    amount = 1.0
    if not isinstance(species, dict):
        water = 0.0
        try:
            masses = elemental_masses(feed)
            water, amount = tidy(masses['water']), tidy(masses['solids'])
        except ValueError:  # FeedError: no usable elemental analysis
            pass
        salts = feed.get('salts')
        feed.clear()
        species = feed['species'] = {'water': {'value': water, 'unit': 'kg'}}
        if salts:
            feed['salts'] = salts
    if name in species:
        raise ValueError(f'{name} is already in the feed.')
    species[name] = {'value': amount, 'unit': 'kg'}
    return 'the feed'


class SpeciesTab:
    def __init__(self, state, editor):
        self.state = state
        self.editor = editor
        self.catalogue = build_catalogue(state.model)
        self.choice = None
        self.dialog_host = ui.element()

        with ui.row().classes('w-full items-start gap-3'):
            self.text = ui.input('Name or formula', on_change=self.fill).props('dense clearable').classes('w-56')
            self.elements = ui.input('Contains elements', placeholder='C, H, O', validation=element_error,
                                     on_change=self.fill).props('dense').classes('w-44')
            self.exact = ui.checkbox('Exact match', on_change=self.fill).classes('mt-2')
            self.source = ui.select(SOURCE_OPTIONS, value='Any', label='Data source', on_change=self.fill) \
                .props('dense').classes('w-36')
            self.show = ui.select(SHOW_OPTIONS, value=SHOW_OPTIONS[0], label='Show', on_change=self.fill) \
                .props('dense').classes('w-40')
            with ui.column().classes('gap-0'):
                ui.label('Can be in').classes('text-caption text-grey-7')
                with ui.row().classes('gap-0'):
                    self.phase_chips = {key: ui.chip(label, selectable=True, color='blue-1', text_color='blue-10',
                                                     on_selection_change=self.fill).props('dense')
                                        for key, label in PHASE_LABELS.items()}
            ui.space()
            self.count = ui.label().classes('text-body2 self-center')

        with ui.row().classes('w-full items-center gap-2'):
            self.selected_label = ui.label('Nothing selected').classes('text-weight-medium')
            ui.label('Add to:').classes('ml-4')
            self.add_buttons = {}
            for target, label in [*PHASE_LABELS.items(), ('feed', 'Feed')]:
                self.add_buttons[target] = ui.button(label, on_click=lambda t=target: self.add(t)) \
                    .props('outline dense no-caps')
            ui.label(ADD_NOTE).classes('text-caption text-grey-7')

        with ui.row().classes('w-full items-center gap-2'):
            self.new_button = ui.button('+ New species', on_click=self.new_species).props('dense no-caps')
            self.edit_button = ui.button('Edit selected', on_click=self.edit_selected) \
                .props('outline dense no-caps')
            self.import_button = ui.button('Import from file…', on_click=self.import_file) \
                .props('outline dense no-caps')
            self.remove_button = ui.button('Remove from my database', on_click=self.confirm_remove) \
                .props('outline dense no-caps color=negative')
            self.database_note = ui.label().classes('text-caption text-grey-7')

        self.table = ui.table(rows=[], columns=COLUMNS, row_key='Species', selection='single', pagination=20,
                              on_select=self.selected).classes('w-full').props('dense flat bordered')
        self.table.add_slot('body-cell-can', CAN_BE_IN_SLOT)
        self.table.add_slot('body-cell-data', DATA_SLOT)
        self.table.add_slot('body-cell-status', STATUS_SLOT)
        self.fill()

    # ------------------------------------------------------------ table
    @property
    def user_dir(self):
        return self.state.model.user_data

    def shown(self):
        """The catalogue rows that pass the filters (unknown element symbols are left out)."""
        symbols = re.split(r'[,;\s]+', self.elements.value or '')
        elements = ', '.join(s for s in symbols if s and not element_error(s))
        return filter_catalogue(self.catalogue, text=self.text.value or '', elements=elements,
                                exact=bool(self.exact.value), source=self.source.value or 'Any',
                                show=self.show.value or 'All',
                                can_be_in=[k for k, chip in self.phase_chips.items() if chip.selected])

    def fill(self, _=None):
        shown = self.shown()
        self.count.set_text(count_text(shown, self.catalogue))
        rows = [{k: row[k] for k in ROW_FIELDS} for row in table_rows(shown)]
        self.table.rows = rows
        self.table.selected = self._rows_for(self.choice['Species'] if self.choice else None)
        self.table.update()
        self.selected()

    def _rows_for(self, name):
        if name is None:
            return []
        return [{k: row[k] for k in ROW_FIELDS}
                for row in table_rows(self.catalogue[self.catalogue.Species == name])]

    def select(self, name):
        """Select ``name`` (also when the filters hide it)."""
        self.table.selected = self._rows_for(name)
        self.table.update()
        self.selected()

    def selected(self, _=None):
        rows = self.table.selected
        self.choice = rows[0] if rows else None
        choice = self.choice
        usable = bool(choice and choice['Usable'])
        phases = {key for key, label in PHASE_LABELS.items() if label in (choice or {}).get('Can be in', [])}
        for target, button in self.add_buttons.items():
            enabled = usable and (choice['Feed'] if target == 'feed' else target in phases)
            button.set_enabled(bool(enabled))
        if choice is None:
            self.selected_label.set_text('Nothing selected')
        else:
            self.selected_label.set_text(f'{choice["Species"]} ({choice["Formula"]})'
                                         if choice['Formula'] else choice['Species'])
        has_db = self.user_dir is not None
        self.database_note.set_text(DATABASE_NOTE if has_db else DATABASE_OFF)
        for button in (self.new_button, self.import_button):
            button.set_enabled(has_db)
        self.edit_button.set_enabled(has_db and choice is not None)
        self.remove_button.set_visibility(bool(has_db and choice and choice['Yours']))

    # ------------------------------------------------------------ add to the case
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

    # ------------------------------------------------------------ the user database
    def new_species(self):
        self.dialog_host.clear()
        with self.dialog_host:
            self.dialog = SpeciesDialog(self.state.model, self.user_dir, self.after_save).open()

    def edit_selected(self):
        if self.choice is None:
            return
        name = self.choice['Species']
        try:
            spec = spec_from_model(self.state.model, name)
        except (KeyError, ValueError) as exc:
            ui.notify(f'{name} cannot be edited here: {exc}', type='warning')
            return
        self.dialog_host.clear()
        with self.dialog_host:
            self.dialog = SpeciesDialog(self.state.model, self.user_dir, self.after_save, spec=spec, editing=name,
                                        overridden=bool(self.choice['Overrides']),
                                        user_only=bool(self.choice['Yours'] and not self.choice['Overrides']))
            self.dialog.open()

    def import_file(self):
        self.dialog_host.clear()
        with self.dialog_host:
            self.dialog = ImportDialog(self.state.model, self.after_import).open()

    def confirm_remove(self):
        if self.choice is None or not self.choice['Yours']:
            return
        name, overrides = self.choice['Species'], self.choice['Overrides']
        self.dialog_host.clear()
        with self.dialog_host, ui.dialog() as dialog, ui.card():
            ui.label(f'Remove your values for {name}? The shipped values are used again.' if overrides
                     else f'Remove {name} from your database?')
            with ui.row().classes('w-full'):
                ui.space()
                ui.button('Cancel', on_click=dialog.close).props('flat no-caps')
                ui.button('Remove', on_click=lambda: self.remove(name, dialog)).props('no-caps color=negative') \
                    .mark('confirm-remove')
        dialog.open()

    async def remove(self, name, dialog=None):
        if dialog is not None:
            dialog.close()
        try:
            UserDatabase.load(self.user_dir).remove_species(name)
        except OSError as exc:
            ui.notify(f'Not removed: {exc}', type='negative')
            return
        await self.apply_changes(f'Removed {name} from your database')

    async def after_save(self, name, action):
        if action == 'reset':
            await self.apply_changes(f'{name} uses the shipped values again')
            return
        await self.apply_changes(None, select=name)
        row = self.catalogue[self.catalogue.Species == name]
        status = row.Status.iloc[0] if len(row) else 'Not listed'
        ui.notify(f'Saved {name}' if status == 'Usable' else f'Saved {name}. Not usable yet: {status}',
                  type='positive' if status == 'Usable' else 'warning', multi_line=True)

    async def after_import(self, count):
        await self.apply_changes(f'Imported {count} species')

    async def apply_changes(self, message, select=None):
        """Use the changed user database: reload the model, the lists and the workers."""
        model = await background.io_bound(Model.load, **self.state.model.load_args)
        self.state.use_model(model)
        self.catalogue = build_catalogue(model)
        now = self.state.pool.restart(model.load_args)
        if select is not None:
            self.choice = {'Species': select}
        self.fill()
        self.editor.form.refresh()
        self.editor.changed()
        if message:
            ui.notify(message, type='positive')
        if not now:
            ui.notify(RESTART_LATER)
