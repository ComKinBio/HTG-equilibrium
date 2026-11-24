"""The Case tab: a form over a case file (CaseDocument.data).

The form writes case-file format 2 only. Everything it computes (point count,
feed numbers, notices, the first problem) comes from pure functions in
``case_form`` and ``htg.feed_v2``, called on every change.
"""
import copy
import math
import re
from pathlib import Path

from nicegui import run as background
from nicegui import ui

from .. import units
from ..constants import ATOMIC_WEIGHTS
from ..feed_v2 import (FeedError, elemental_masses, feed_mode, molar_mass_function, registry_elements, scaled_elements,
                       water_molar_mass)
from .case_form import (ADVANCED_CAPTION, AQUEOUS_MODELS, ATOMIC_WEIGHTS_CAPTION, CARBON_CAPTION, CARBON_CHOICES,
                        CARBON_PHASES, CARBON_RANGE, CARBON_RANGE_LABEL, CHECK_CAPTION, COMPOSITION_TOL_CAPTION,
                        CONTINUATION_TOOLTIP, CUTOFF_CAPTION, DISTINCT_CAPTION, ELEMENTS_CAPTION, ENERGY_SCALE_CAPTION,
                        FD_METHOD, FEED_CAPTION, FEED_CAPTION_SPECIES, FEED_TOGGLE, GIBBS_TOL_CAPTION,
                        INITIAL_COMPOSITION_CAPTION, INITIAL_COMPOSITIONS, INPUT_FORMATS, MEASURED_CAPTION,
                        METHOD_TOOLTIP, METHODS, MODE_TOGGLE,
                        PHASE_NAMES, PRODUCTS_CAPTION, QUANTITY_CHOICES, SEED_CAPTION, SOLID_CARBON, SPECIES_FOOTER,
                        STABILITY_TEST, carbon_bound, carbon_choice, check_feed, conflict_notice, convert_share,
                        convert_species_amount, element_error, element_symbol, entered_sum_notice, file_status,
                        format_kg, format_mol, format_percent, get, initial_guesses_tooltip, number_text, parse_number,
                        point_count, point_count_text, problem_text, put, ready_text, tidy)
from .document import LEGACY_NOTICE, SALTS_NOTICE, CaseDocument, clean

__all__ = ['CaseEditor', 'PHASE_NAMES', 'QUANTITY_CHOICES', 'get', 'put', 'parse_numbers', 'parse_points',
           'format_numbers']

MAX_LIST_VALUES = 10000
_NUMBER_FORMATS = 'Write numbers like 300, 350, a range 300:400:50 or 300:400 n=5.'
# order of the form's sections, for "the first problem" in the status line
SECTIONS = ('name', 'conditions', 'feed', 'products', 'carbon', 'aqueous', 'solver', 'advanced')
TEMPERATURE_UNITS = {u: units.label(u) for u in units.TEMPERATURE_UNITS}
PRESSURE_UNITS = {u: units.label(u) for u in units.PRESSURE_UNITS}
SHARE_UNITS = {'kg': 'kg', '%': '%'}
SPECIES_UNITS = {'mol': 'mol', 'kg': 'kg'}
CAPTION = 'text-caption text-grey-7'
TITLE = 'text-subtitle1 text-weight-medium'


# ---------------------------------------------------------------- number lists
def _number(text):
    try:
        value = float(text)
    except ValueError:
        value = math.nan
    if not math.isfinite(value):
        raise ValueError(f'"{text}" is not a number. {_NUMBER_FORMATS}')
    return value


def _clean(value):
    """Drop floating-point noise from generated or converted values: 299.99999999999994 -> 300."""
    return float(f'{value:.12g}')


def _range(token, count):
    """'a:b:step' (count None) or 'a:b' with count -> values (in the unit typed)."""
    parts = token.split(':')
    if len(parts) not in (2, 3) or '' in parts:
        raise ValueError(f'"{token}" is not a range. Write start:stop:step (300:400:50) '
                         'or start:stop n=count (300:400 n=5).')
    start, stop = _number(parts[0]), _number(parts[1])
    if start > stop:
        raise ValueError(f'In "{token}" the start is above the stop. Put the smaller value first.')
    if len(parts) == 3:
        if count is not None:
            raise ValueError(f'"{token} n={count}" gives both a step and a count. Use one of them.')
        step = _number(parts[2])
        if step <= 0:
            raise ValueError(f'In "{token}" the step must be above zero.')
        steps = int(round((stop - start) / step)) if stop > start else 0
        if steps + 1 > MAX_LIST_VALUES:
            raise ValueError(f'"{token}" gives more than {MAX_LIST_VALUES} values. Use a larger step.')
        return [_clean(start + i * step) for i in range(steps + 1) if start + i * step <= stop + 1e-9 * step]
    if count is None:
        raise ValueError(f'"{token}" needs a step or a count, e.g. {token}:10 or {token} n=5.')
    return [_clean(start + i * (stop - start) / (count - 1)) for i in range(count - 1)] + [stop]


def _count(token):
    text = token[2:]
    if not text.isdigit() or int(text) < 2:
        raise ValueError(f'"{token}": the count after n= must be a whole number of 2 or more.')
    if int(text) > MAX_LIST_VALUES:
        raise ValueError(f'"{token}": the count can be at most {MAX_LIST_VALUES}.')
    return int(text)


def parse_numbers(text, unit=None, sort=True):
    """Text of a number field -> list of floats.

    Accepted (separators: comma, semicolon or whitespace; formats can be mixed):

        300, 350, 400     list                          300, 350, 400
        300:400:50        range start:stop:step         300, 350, 400 (stops at or before stop)
        300:400 n=5       evenly spaced start:stop n=count   300, 325, 350, 375, 400
        250, 300:400:50   mixed                         250, 300, 350, 400

    ``n=`` belongs to the ``a:b`` range just before it; the count is a whole
    number of 2 or more.

    unit  the unit the numbers are typed in (a key of ``htg.units``, e.g. 'K' or
          'MPa'); the result is converted to the file unit (°C or bar). None:
          the numbers are already in file units.
    sort  sort and remove duplicates (temperature and pressure grids). False
          keeps the order as typed (e.g. one 'T, P' point).

    Raises ValueError with a message that can be shown under the field.
    """
    text = re.sub(r'\s*:\s*', ':', (text or '').strip())
    text = re.sub(r'\bn\s*=\s*', 'n=', text)
    if not text:
        return []
    tokens = [t for t in re.split(r'[,;\s]+', text) if t]
    values = []
    for i, token in enumerate(tokens):
        if token.startswith('n='):
            previous = tokens[i - 1] if i else ''
            if previous.count(':') != 1:
                raise ValueError(f'"{token}" must follow a range start:stop, e.g. 300:400 n=5.')
            continue
        if ':' in token:
            following = tokens[i + 1] if i + 1 < len(tokens) else ''
            count = _count(following) if following.startswith('n=') else None
            values += _range(token, count)
        else:
            values.append(_number(token))
    if unit is not None and units.unit_key(unit) not in units.FILE_UNITS.values():
        values = [_clean(units.to_file(v, unit)) for v in values]
    return sorted(set(values)) if sort else values


def format_numbers(value, unit=None):
    """Stored value(s) in °C or bar -> field text in ``unit`` (None: as stored)."""
    if value is None:
        return ''
    values = value if isinstance(value, list) else [value]
    if unit is None:
        return units.format_list(values)
    return ', '.join(units.display(v, unit) for v in values)


def parse_points(text, temperature_unit=None, pressure_unit=None):
    """Lines 'T, P' (in the given units) -> [[T_C, P_bar], ...]."""
    points = []
    for line in (text or '').splitlines():
        if line.strip():
            numbers = parse_numbers(line, sort=False)
            if len(numbers) != 2:
                raise ValueError(f'"{line.strip()}": give a temperature and a pressure.')
            T, P = numbers
            if temperature_unit is not None:
                T = _clean(units.to_file(T, temperature_unit))
            if pressure_unit is not None:
                P = _clean(units.to_file(P, pressure_unit))
            points.append([T, P])
    return points


def scalar(value):
    return value[0] if isinstance(value, list) and len(value) == 1 else value


def default_feed(mode, feed_species=()):
    """A starting feed for a mode (the template's values)."""
    if mode == 'species':
        species = {'water': {'value': 4.5, 'unit': 'kg'}}
        if 'cellulose' in feed_species:
            species['cellulose'] = {'value': 0.5, 'unit': 'kg'}
        return {'species': species}
    return {'total_kg': 5.0, 'solids': {'value': 20.0, 'unit': '%'}, 'ash': {'value': 8.44, 'unit': '%'},
            'elements_wt_percent': {'C': 30.4, 'H': 4.41, 'O': 53.58, 'N': 2.66, 'S': 0.51}}


# ---------------------------------------------------------------- editor
class CaseEditor:
    """File bar, status line and the form of the Case tab.

    Text fields are parsed on every change. A field that cannot be parsed keeps
    its text, shows the message under it and leaves the case data unchanged;
    any problem disables Save and Add to queue and is named in the status line.
    Old files that cannot be converted (legacy_blocker) are not shown in the
    form: they keep the old format and can only be checked, queued or saved as
    they are.
    """

    def __init__(self, state, on_queued):
        self.state = state
        self.on_queued = on_queued
        self.settings = state.settings
        self.fields = {}        # key -> (element, label, parse, store)
        self.parse_errors = {}  # key -> message
        self.feed_errors = {}   # key -> message from the feed check
        self.feed_stash = {}    # feed of the other mode, to switch back
        self.built = False
        self._order = 0
        self.snapshot = None
        with ui.row().classes('w-full items-center gap-2 no-wrap'):
            self.file_select = ui.select(self._case_files(), label='Case file', with_input=True,
                                         on_change=self._open_selected).classes('w-72')
            ui.button('New', icon='note_add', on_click=self._new).props('flat no-caps')
            self.save_button = ui.button('Save', icon='save', on_click=self._save).props('flat no-caps') \
                .mark('case-save')
            self.save_as_button = ui.button('Save as…', icon='save_as', on_click=self._save_as) \
                .props('flat no-caps').mark('case-save-as')
            self.file_label = ui.label().classes(CAPTION)
            ui.space()
            ui.button('Check', icon='fact_check', on_click=self._check).props('outline no-caps')
            self.queue_button = ui.button('Add to queue', icon='play_arrow', on_click=self._queue) \
                .props('color=primary no-caps').mark('case-queue')
        self.status = ui.label().classes('text-body2').mark('case-status')
        self.notices = Part(self._build_notices)
        self.form = Part(self._build_form)
        self._take_snapshot()
        self.changed()

    # ------------------------------------------------------------ files
    @property
    def document(self):
        return self.state.document

    @property
    def data(self):
        return self.state.document.data

    @property
    def t_unit(self):
        return self.settings.get('temperature_unit')

    @property
    def p_unit(self):
        return self.settings.get('pressure_unit')

    def _case_files(self):
        root = self.state.cases_dir
        files = sorted(root.glob('**/*.toml')) if root.exists() else []
        return {str(p): str(p.relative_to(root)) for p in files}

    def _molar_mass(self):
        return molar_mass_function(registry_elements(self.state.model.registry))

    def _open_selected(self, event):
        if event.value:
            self.open(event.value)

    def open(self, path):
        try:
            self.state.document = CaseDocument.open(path, molar_mass=self._molar_mass())
        except Exception as exc:
            ui.notify(f'Could not open {path}: {exc}', type='negative', multi_line=True)
            return
        self._reload()

    def _new(self):
        self.state.document = CaseDocument.new()
        self.file_select.set_value(None)
        self._reload()

    def _reload(self):
        self.feed_stash = {}
        self._take_snapshot()
        self.notices.refresh()
        self.form.refresh()
        self.changed()

    def _take_snapshot(self):
        self.snapshot = copy.deepcopy(clean(self.data))

    def dirty(self):
        return clean(self.data) != self.snapshot

    def _save(self):
        if self.document.path is None:
            return self._save_as()
        self._write(self.document.path)

    def _save_as(self):
        default = self.state.cases_dir / f"{self.data.get('name') or 'case'}.toml"
        with ui.dialog() as dialog, ui.card().classes('w-[36rem]'):
            ui.label('Save case file as').classes('text-h6')
            path = ui.input('File', value=str(default)).classes('w-full')
            with ui.row():
                ui.button('Save', on_click=lambda: (self._write(Path(path.value)), dialog.close()))
                ui.button('Cancel', on_click=dialog.close).props('flat')
        dialog.open()

    def _write(self, path):
        problem = self.first_problem()
        if problem:
            ui.notify(f'Not saved: {problem}', type='negative', multi_line=True)
            return
        try:
            saved = self.document.save(path)
        except Exception as exc:
            ui.notify(f'Not saved: {exc}', type='negative', multi_line=True)
            return
        ui.notify(f'Saved {saved}', type='positive')
        self.file_select.set_options(self._case_files())
        self._take_snapshot()
        self.notices.refresh()
        self.changed()

    def _build_notices(self):
        document = self.document
        for text in ([document.legacy_blocker] if document.legacy_blocker else []) + \
                ([LEGACY_NOTICE] if document.converted_from_legacy else []):
            with ui.row().classes('w-full items-center no-wrap bg-amber-1 q-pa-sm rounded-borders'):
                ui.icon('info', color='warning')
                ui.label(text).classes('text-body2').mark('case-notice')

    # ------------------------------------------------------------ status and actions
    def set(self, path, value):
        put(self.data, path, value)
        self.changed()

    def problems(self):
        """Problems with the form, first one first (field problems in form order, then the case)."""
        out = []
        if self.built:
            for key, (_, label, _, _) in sorted(self.fields.items(), key=lambda item: self._rank(item[0])):
                message = self.parse_errors.get(key) or self.feed_errors.get(key)
                if message:
                    out.append(problem_text(label, message))
            if self.feed_message and not any(self.feed_message == v for v in self.feed_errors.values()):
                out.append(self.feed_message)
        return out

    def first_problem(self):
        problems = self.problems()
        if problems:
            return problems[0]
        try:
            self.document.case()
        except Exception as exc:
            return str(exc)
        return None

    def _rank(self, key):
        section = key.split('.')[0]
        return SECTIONS.index(section) if section in SECTIONS else len(SECTIONS), self.fields[key][0].id

    def changed(self):
        """Update everything derived from the data: numbers, notices, field messages, status."""
        if self.built:
            self._update_derived()
        self.file_label.set_text(file_status(self.document.path, self.dirty()))
        problems = self.problems()
        case = None
        if not problems:
            try:
                case = self.document.case()
            except Exception as exc:
                problems = [str(exc)]
        if problems:
            self.status.set_text(f'✗ {problems[0]}')
            self.status.classes(replace='text-body2 text-negative')
        else:
            self.status.set_text(ready_text(len(case.grid()), len(case.chains())))
            self.status.classes(replace='text-body2 text-positive')
        for button in (self.save_button, self.save_as_button, self.queue_button):
            button.set_enabled(not problems)

    async def _check(self):
        problem = self.first_problem()
        if problem:
            ui.notify(problem, type='negative', multi_line=True)
            return
        from ..sweep import preflight_messages
        case = self.document.case()
        notification = ui.notification('Checking every point…', spinner=True, timeout=None)
        try:
            messages = await background.io_bound(preflight_messages, self.state.model, case)
        finally:
            notification.dismiss()
        with ui.dialog() as dialog, ui.card().classes('w-[48rem] max-w-full'):
            if messages:
                ui.label('Problems found').classes('text-negative text-h6')
                with ui.column().classes('gap-1'):
                    for message in messages[:50]:
                        ui.label(message).classes('text-body2')
                    if len(messages) > 50:
                        ui.label(f'… and {len(messages) - 50} more.').classes('text-body2')
            else:
                ui.label(f'All {point_count_text(len(case.grid()))} passed the checks.') \
                    .classes('text-positive text-h6')
            ui.label(CHECK_CAPTION).classes('text-caption')
            ui.button('Close', on_click=dialog.close).props('flat')
        dialog.open()

    def _queue(self):
        problem = self.first_problem()
        if problem:
            ui.notify(f'Not queued: {problem}', type='negative', multi_line=True)
            return
        try:
            case = self.document.case()
            run = self.state.submit(case)
        except Exception as exc:
            ui.notify(f'Not queued: {exc}', type='negative', multi_line=True)
            return
        ui.notify(f'Queued {case.name}: {point_count_text(run.total)}', type='positive')
        self.on_queued(run)

    # ------------------------------------------------------------ fields
    def _field(self, key, label, text, parse, store, *, field_label=None, textarea=False, classes='w-56',
               placeholder=None, hint=None):
        """A text field parsed on every change: parse(text) -> value or ValueError (shown
        under the field), then store(value) writes the case data."""
        maker = ui.textarea if textarea else ui.input
        shown = label if field_label is None else field_label
        element = maker(shown or None, value=text).classes(classes).mark(f'case-{key}').without_auto_validation()
        if placeholder is not None:
            element.props(f'placeholder="{placeholder}" stack-label')
        if hint:
            set_hint(element, hint)
        if textarea:
            element.props('autogrow')
        if field_label == '':  # a cell of a table
            element.props('dense hide-bottom-space')
        self.fields[key] = (element, label, parse, store)
        element.on_value_change(lambda _: self.apply(key))
        return element

    def apply(self, key):
        """Parse the field's text and store it (or record the problem), then update the form."""
        element, _, parse, store = self.fields[key]
        try:
            value = parse(element.value)
        except ValueError as exc:
            self.parse_errors[key] = str(exc)
        else:
            self.parse_errors.pop(key, None)
            store(value)
        self.changed()

    def _forget(self, prefix):
        """Drop the fields (and their problems) of a part of the form that is rebuilt."""
        for key in [k for k in self.fields if k.startswith(prefix)]:
            del self.fields[key]
            self.parse_errors.pop(key, None)

    def _optional(self, key, label, path, default, caption=None, **limits):
        """A number field that may be left empty (= the default, shown as placeholder)."""
        def store(value):
            put(self.data, path, value)
        element = self._field(key, label, number_text(get(self.data, path)),
                              lambda text: parse_number(text, **limits), store, classes='w-72',
                              placeholder=default)
        if caption:
            ui.label(caption).classes(CAPTION + ' w-72')
        return element

    # ------------------------------------------------------------ the form
    def _build_form(self):
        self.fields = {}
        self.parse_errors = {}
        self.feed_errors = {}
        self.feed_message = None
        self.built = False
        if self.document.legacy_blocker:
            return
        data = self.data
        with ui.column().classes('w-full gap-3'):
            with ui.card().classes('w-full'):
                with ui.row().classes('w-full no-wrap items-start'):
                    self._field('name.name', 'Name (also the results folder name)', data.get('name', ''),
                                lambda text: text.strip(), lambda value: put(self.data, ['name'], value),
                                classes='w-80')
                    ui.textarea('Notes', value=data.get('description', ''),
                                on_change=lambda e: self.set(['description'], e.value)) \
                        .classes('grow').props('autogrow rows=1')
            with card('Temperature and pressure'):
                self._conditions = Part(self._build_conditions)
            if self.document.has_salts:
                with ui.row().classes('w-full items-center no-wrap bg-amber-1 q-pa-sm rounded-borders'):
                    ui.icon('info', color='warning')
                    ui.label(SALTS_NOTICE).classes('text-body2').mark('case-notice')
            with card('Feed'):
                self._feed = Part(self._build_feed)
            with card('Possible products', PRODUCTS_CAPTION):
                self._products()
            with ui.row().classes('w-full items-stretch no-wrap gap-3'):
                with card('Carbon split', CARBON_CAPTION, grow=True):
                    self._carbon()
                with card('Aqueous phase', grow=True):
                    self._aqueous()
                with card('Solver', grow=True):
                    self._solver()
            self._advanced()
        self.built = True

    # ------------------------------------------------------------ temperature and pressure
    def _build_conditions(self):
        self._forget('conditions.')
        T, P = self.t_unit, self.p_unit
        conditions = self.data.setdefault('conditions', {})
        explicit = conditions.get('points') is not None
        with ui.row().classes('w-full items-center gap-4'):
            ui.toggle(MODE_TOGGLE, value='points' if explicit else 'grid',
                      on_change=lambda e: self._switch_conditions(e.value)).props('dense no-caps') \
                .mark('case-conditions-mode')
            ui.select(TEMPERATURE_UNITS, label='Temperature unit', value=T,
                      on_change=lambda e: self._set_unit('temperature_unit', e.value)).classes('w-36') \
                .mark('case-temperature-unit')
            ui.select(PRESSURE_UNITS, label='Pressure unit', value=P,
                      on_change=lambda e: self._set_unit('pressure_unit', e.value)).classes('w-36') \
                .mark('case-pressure-unit')
        with ui.row().classes('w-full items-start gap-4 no-wrap'):
            if explicit:
                text = '\n'.join(f'{units.display(t, T)}, {units.display(p, P)}' for t, p in conditions['points'])
                self._field('conditions.points', f'(T, P) points, one per line [{units.label(T)}, {units.label(P)}]',
                            text, self._parse_points, lambda v: put(self.data, ['conditions', 'points'], v),
                            textarea=True, classes='grow')
            else:
                self._field('conditions.temperatures', units.field_label('Temperatures', T),
                            format_numbers(conditions.get('temperature_C'), T),
                            lambda text: self._parse_list(text, self.t_unit),
                            lambda v: put(self.data, ['conditions', 'temperature_C'], scalar(v)), classes='grow')
                self._field('conditions.pressures', units.field_label('Pressures', P),
                            format_numbers(conditions.get('pressure_bar'), P),
                            lambda text: self._parse_list(text, self.p_unit),
                            lambda v: put(self.data, ['conditions', 'pressure_bar'], scalar(v)), classes='w-80')
            self.point_label = ui.label().classes('text-body2 q-mt-md text-no-wrap').mark('case-point-count')
        with ui.expansion('Input formats', value=self.settings.get('input_formats_open'),
                          on_value_change=lambda e: self.settings.set('input_formats_open', bool(e.value))) \
                .props('dense').classes('w-full').mark('case-input-formats'):
            ui.table(columns=[{'name': c.lower(), 'label': c, 'field': c.lower(), 'align': 'left'}
                              for c in ('Format', 'Type', 'Gives')],
                     rows=[{'format': f, 'type': t, 'gives': g} for f, t, g in INPUT_FORMATS],
                     row_key='format').props('dense flat')

    @staticmethod
    def _parse_list(text, unit):
        values = parse_numbers(text, unit)
        if not values:
            raise ValueError('Give at least one value.')
        return values

    def _parse_points(self, text):
        points = parse_points(text, self.t_unit, self.p_unit)
        if not points:
            raise ValueError('Give at least one point.')
        return points

    def _set_unit(self, key, unit):
        if unit is None or unit == self.settings.get(key):
            return
        self.settings.set(key, unit)
        self._conditions.refresh()
        if key == 'temperature_unit' and getattr(self, 'cutoff_field', None) is not None:
            field = self.cutoff_field
            field._props['label'] = units.field_label('Remove at and above', unit)
            if 'aqueous.cutoff' not in self.parse_errors:
                field.set_value(format_numbers(get(self.data, ['aqueous', 'cutoff_C']), unit))
            field.update()
        self.changed()

    def _switch_conditions(self, mode):
        conditions = self.data.setdefault('conditions', {})
        if mode == 'points' and 'points' not in conditions:
            temperatures = conditions.get('temperature_C')
            pressures = conditions.get('pressure_bar')
            temperatures = temperatures if isinstance(temperatures, list) else [temperatures]
            pressures = pressures if isinstance(pressures, list) else [pressures]
            points = [[t, p] for p in pressures for t in temperatures if t is not None and p is not None]
            conditions.clear()
            conditions['points'] = points or [[300.0, 240.0]]
        elif mode == 'grid' and 'points' in conditions:
            points = conditions.pop('points') or [[300.0, 240.0]]
            conditions['temperature_C'] = sorted({t for t, _ in points})
            conditions['pressure_bar'] = scalar(sorted({p for _, p in points}))
        self._conditions.refresh()
        self.changed()

    # ------------------------------------------------------------ feed
    def _feed_data(self):
        return self.data.setdefault('feed', {})

    def _build_feed(self):
        self._forget('feed.')
        feed = self._feed_data()
        mode = feed_mode(feed)
        self.feed_mode = mode
        self.used_labels = {}
        ui.label(FEED_CAPTION_SPECIES if mode == 'species' else FEED_CAPTION).classes(CAPTION)
        ui.toggle(FEED_TOGGLE, value=mode, on_change=lambda e: self._switch_feed(e.value)) \
            .props('dense no-caps').mark('case-feed-mode')
        with ui.row().classes('w-full items-start gap-10 no-wrap'):
            with ui.column().classes('gap-2'):
                if mode == 'species':
                    self._species_feed()
                else:
                    self._elemental_feed()
            with ui.column().classes('gap-1 q-ml-auto'):
                ui.label('Feed summary').classes('text-weight-medium')
                self._summary = Part(self._build_summary, 'gap-0')

    def _switch_feed(self, mode):
        feed = self._feed_data()
        current = feed_mode(feed)
        if mode is None or mode == current:
            return
        self.feed_stash[current] = {k: copy.deepcopy(v) for k, v in feed.items() if k != 'salts'}
        new = copy.deepcopy(self.feed_stash.get(mode) or default_feed(mode, self.state.feed_species))
        if feed.get('salts'):
            new['salts'] = feed['salts']
        self.data['feed'] = new
        self._feed.refresh()
        self.changed()

    def _masses(self):
        try:
            return elemental_masses(self._feed_data())
        except (FeedError, ValueError):
            return None

    def _elemental_feed(self):
        feed = self._feed_data()
        self._field('feed.total', '1 · Total mass [kg]', number_text(feed.get('total_kg')),
                    lambda text: parse_number(text, required=True, above=0),
                    lambda v: put(self._feed_data(), ['total_kg'], v), hint='Solids + water')
        self._share('solids', '2 · Solids', required=True)
        self._share('ash', '3 · Ash', required=False)
        self._elements = Part(self._build_elements)

    def _share(self, name, label, required):
        """Solids or ash: a value with a kg | % switch."""
        entry = self._feed_data().get(name) or {}
        unit = {'unit': entry.get('unit', '%')}

        def hint():
            if name == 'solids':
                return '% of total mass' if unit['unit'] == '%' else 'kg'
            return ('% of solids' if unit['unit'] == '%' else 'kg') + '. Empty = no ash.'

        def store(value):
            put(self._feed_data(), [name], None if value is None else {'value': value, 'unit': unit['unit']})

        def switch(e):
            old, unit['unit'] = unit['unit'], e.value
            set_hint(element, hint())
            masses = self._masses()
            whole = None if masses is None else masses['total'] if name == 'solids' else masses['solids']
            current = (self._feed_data().get(name) or {}).get('value')
            if f'feed.{name}' in self.parse_errors or current is None:
                self.apply(f'feed.{name}')
                return
            converted = convert_share(current, old, unit['unit'], whole)
            if converted == current:
                self.apply(f'feed.{name}')
            else:
                element.set_value(number_text(converted))  # stores in the new unit

        with ui.row().classes('items-start no-wrap gap-2'):
            element = self._field(f'feed.{name}', label, number_text(entry.get('value')),
                                  lambda text: parse_number(text, required=required,
                                                            **({'above': 0} if required else {'minimum': 0})),
                                  store, hint=hint())
            ui.toggle(SHARE_UNITS, value=unit['unit'], on_change=switch).props('dense no-caps') \
                .classes('q-mt-sm').mark(f'case-feed.{name}-unit')

    def _build_elements(self):
        self._forget('feed.element:')
        elements = self._feed_data().get('elements_wt_percent') or {}
        self.used_labels = {}
        with ui.row().classes('w-full items-baseline no-wrap gap-6 q-mt-sm'):
            ui.label('4 · Elements [mass %, ash-free solids]')
            ui.space()
            ui.label(ELEMENTS_CAPTION).classes(CAPTION)
        with ui.grid(columns='3rem 7rem 5rem 2rem').classes('items-center gap-x-4 gap-y-0'):
            for text in ('Element', 'Entered', 'Used', ''):
                ui.label(text).classes(CAPTION)
            for element in list(elements):
                ui.label(element)
                self._field(f'feed.element:{element}', f'Element {element}', number_text(elements[element]),
                            lambda text: parse_number(text, required=True, minimum=0),
                            lambda v, e=element: self._feed_data().setdefault('elements_wt_percent', {})
                            .__setitem__(e, v), field_label='', classes='w-28')
                self.used_labels[element] = ui.label().mark(f'case-used-{element}')
                ui.button(icon='close', on_click=lambda e=element: self._remove_element(e)) \
                    .props('flat dense round size=sm').mark(f'case-remove-{element}')
            ui.label('Sum').classes('text-weight-medium')
            self.entered_sum_label = ui.label().classes('text-weight-medium').mark('case-entered-sum')
            self.used_sum_label = ui.label().classes('text-weight-medium')
            ui.label('')
        with ui.row().classes('items-start gap-2'):
            symbol = ui.input('Element symbol').classes('w-32').mark('case-new-element')
            ui.button('+ Add element', on_click=lambda: self._add_element(symbol)) \
                .props('flat dense no-caps').classes('q-mt-sm').mark('case-add-element')
        self.scale_notice = ui.label().classes('text-body2 text-warning').mark('case-scale-notice')
        self.feed_general = ui.label().classes('text-caption text-negative')

    def _add_element(self, symbol_input):
        elements = self._feed_data().setdefault('elements_wt_percent', {})
        symbol = element_symbol(symbol_input.value)
        error = element_error(symbol, elements)
        symbol_input.error = error
        if error:
            return
        elements[symbol] = 0.0
        self._elements.refresh()
        self.changed()

    def _remove_element(self, element):
        (self._feed_data().get('elements_wt_percent') or {}).pop(element, None)
        self._elements.refresh()
        self.changed()

    # species amounts
    def _species_elements(self):
        return registry_elements(self.state.model.registry)

    def _species_molar_mass(self, name):
        try:
            return self._molar_mass()(name)
        except (FeedError, ValueError, KeyError):
            return None

    def _species_feed(self):
        species = self._feed_data().setdefault('species', {})
        if 'water' not in species:
            self._feed_data()['species'] = species = {'water': {'value': 0.0, 'unit': 'kg'}, **species}
        self.species_items = [[name, entry.get('value'), entry.get('unit', 'kg')]
                              for name, entry in species.items() if name != 'water']
        water = species['water']
        water_unit = {'unit': water.get('unit', 'kg')}

        def store_water(value):
            self._feed_data()['species']['water'] = {'value': value, 'unit': water_unit['unit']}

        def switch_water(e):
            old, water_unit['unit'] = water_unit['unit'], e.value
            current = self._feed_data()['species']['water'].get('value')
            converted = convert_species_amount(current, old, e.value,
                                               water_molar_mass(self.data.get('atomic_weights')))
            if 'feed.water' in self.parse_errors or converted == current:
                self.apply('feed.water')
            else:
                water_field.set_value(number_text(converted))

        with ui.row().classes('items-start no-wrap gap-2'):
            ui.label('Water (always present)').classes('w-64 q-mt-md')
            water_field = self._field('feed.water', 'Water', number_text(water.get('value')),
                                      lambda text: parse_number(text, required=True, minimum=0), store_water,
                                      field_label='Amount', classes='w-28')
            ui.toggle(SPECIES_UNITS, value=water_unit['unit'], on_change=switch_water) \
                .props('dense no-caps').classes('q-mt-sm').mark('case-feed.water-unit')
        self._species_rows = Part(self._build_species_rows)
        ui.button('+ Add species', on_click=self._add_species).props('flat dense no-caps') \
            .mark('case-add-species')
        self.feed_general = ui.label().classes('text-caption text-negative')

    def _build_species_rows(self):
        self._forget('feed.species')
        for i, (name, value, unit) in enumerate(self.species_items):
            with ui.row().classes('items-start no-wrap gap-2'):
                options = sorted(set(self.state.feed_species) | ({name} if name else set()))
                select = ui.select(options, value=name or None, with_input=True, label='Species',
                                   on_change=lambda e, i=i: self._species_name(i, e.value)) \
                    .classes('w-64').mark(f'case-species-{i}')
                with select.add_slot('after-options'):
                    ui.label(SPECIES_FOOTER).classes(CAPTION + ' q-pa-sm')
                self.fields[f'feed.species-name:{i}'] = (select.without_auto_validation(), 'Species', None, None)
                element = self._field(f'feed.species:{i}', f'Species {name or i + 1}', number_text(value),
                                      lambda text: parse_number(text, required=True, minimum=0),
                                      lambda v, i=i: self._species_value(i, v), field_label='Amount',
                                      classes='w-28')
                ui.toggle(SPECIES_UNITS, value=unit,
                          on_change=lambda e, i=i, f=element: self._species_unit(i, e.value, f)) \
                    .props('dense no-caps').classes('q-mt-sm')
                ui.button(icon='close', on_click=lambda i=i: self._remove_species(i)) \
                    .props('flat dense round size=sm').classes('q-mt-sm')
        self._commit_species()

    def _species_name(self, i, name):
        self.species_items[i][0] = name or ''
        self._commit_species()
        self.changed()

    def _species_value(self, i, value):
        self.species_items[i][1] = value
        self._commit_species()

    def _species_unit(self, i, unit, element):
        name, value, old = self.species_items[i]
        self.species_items[i][2] = unit
        converted = convert_species_amount(value, old, unit, self._species_molar_mass(name) if name else None)
        if f'feed.species:{i}' in self.parse_errors or converted == value:
            self.apply(f'feed.species:{i}')
        else:
            element.set_value(number_text(converted))

    def _add_species(self):
        self.species_items.append(['', 1.0, 'kg'])
        self._species_rows.refresh()
        self.changed()

    def _remove_species(self, i):
        self.species_items.pop(i)
        self._species_rows.refresh()
        self.changed()

    def _commit_species(self):
        """Write the species rows to the feed; unnamed or repeated species are problems of their row."""
        species = self._feed_data().setdefault('species', {})
        out = {'water': species.get('water', {'value': 0.0, 'unit': 'kg'})}
        for i, (name, value, unit) in enumerate(self.species_items):
            key = f'feed.species-name:{i}'
            if not name:
                self.parse_errors[key] = 'Choose a species or remove the row.'
            elif name in out:
                self.parse_errors[key] = f'{name} is already in the feed.'
            else:
                self.parse_errors.pop(key, None)
                if value is not None:
                    out[name] = {'value': value, 'unit': unit}
        self._feed_data()['species'] = out

    # summary
    def _build_summary(self):
        summary = getattr(self, 'feed_summary', None)
        if summary is None:
            ui.label(self.feed_message or '').classes(CAPTION + ' w-64')
            return
        with ui.grid(columns='auto auto auto').classes('gap-x-6 gap-y-0').mark('case-feed-summary'):
            for text in ('', 'kg', 'mol'):
                ui.label(text).classes(CAPTION + ' text-right')
            for row in summary['rows']:
                ui.label(row['name'])
                ui.label(format_kg(row['kg'])).classes('text-right')
                ui.label(format_mol(row['mol'])).classes('text-right')

    def _update_feed(self):
        data = self.data
        feed = self._feed_data()
        species_elements = self._species_elements() if feed_mode(feed) == 'species' else None
        self.feed_summary, field, message = check_feed(feed, data.get('atomic_weights'), species_elements)
        self.feed_errors = {}
        self.feed_message = message
        if message:
            key = {'total': 'feed.total', 'solids': 'feed.solids', 'ash': 'feed.ash',
                   'species:water': 'feed.water'}.get(field)
            if field.startswith('element:'):
                key = 'feed.' + field
            elif field.startswith('species:') and key is None:
                name = field.split(':', 1)[1]
                index = next((i for i, item in enumerate(getattr(self, 'species_items', [])) if item[0] == name), None)
                key = None if index is None else f'feed.species:{index}'
            if key in self.fields:
                self.feed_errors[key] = message
        if self.feed_mode == 'elements':
            try:
                used, entered = scaled_elements(feed.get('elements_wt_percent') or {})
            except FeedError:
                used, entered = {}, None
            for element, label in self.used_labels.items():
                label.set_text(format_percent(used[element]) if element in used else '—')
            self.entered_sum_label.set_text('—' if entered is None else format_percent(entered))
            self.used_sum_label.set_text(format_percent(100.0) if used else '—')
            notice = None if entered is None else entered_sum_notice(entered)
            self.scale_notice.set_text(notice or '')
            self.scale_notice.set_visibility(bool(notice))
        general = message if message and not self.feed_errors else ''
        self.feed_general.set_text(general)
        self.feed_general.set_visibility(bool(general))
        self._summary.refresh()

    # ------------------------------------------------------------ products
    def _products(self):
        for phase, label in PHASE_NAMES.items():
            current = list(get(self.data, ['phases', phase]) or [])
            options = list(dict.fromkeys(self.state.phase_options[phase] + current))
            ui.select(options, label=label, value=current, multiple=True, with_input=True,
                      on_change=lambda e, p=phase: self.set(['phases', p], list(e.value or []))) \
                .classes('w-full').props('use-chips').mark(f'case-phase-{phase}')

    # ------------------------------------------------------------ carbon split
    def _carbon(self):
        for phase in CARBON_PHASES:
            with ui.row().classes('items-start no-wrap gap-2'):
                self._carbon_row(phase)
        target = get(self.data, ['constraints', 'solid_carbon_target'], 'char')
        ui.select(SOLID_CARBON, label='Solid carbon means', value=target,
                  on_change=lambda e: self.set(['constraints', 'solid_carbon_target'],
                                               None if e.value == 'char' else e.value)).classes('w-64')

    def _carbon_row(self, phase):
        choice, first, second = carbon_choice(get(self.data, ['constraints', 'carbon', phase]))
        options = dict(CARBON_CHOICES)
        if choice == CARBON_RANGE:
            options[CARBON_RANGE] = CARBON_RANGE_LABEL
        values = {'first': first, 'second': second}
        name = PHASE_NAMES[phase]
        ui.label(name).classes('q-mt-md w-28')
        select = ui.select(options, value=choice).classes('w-28').mark(f'case-carbon-{phase}')

        def commit():
            kind = select.value
            needed = [values['first']] + ([values['second']] if kind == CARBON_RANGE else [])
            bound = None if kind == 'free' or None in needed else carbon_bound(kind, values['first'], values['second'])
            put(self.data, ['constraints', 'carbon', phase], bound)

        def parser(slot):
            def parse(text):
                if select.value == 'free' or (slot == 'second' and select.value != CARBON_RANGE):
                    return None
                return parse_number(text, required=True, minimum=0, maximum=100)
            return parse

        def storer(slot):
            def store(value):
                values[slot] = None if value is None else tidy(units.percent_to_fraction(value))
                commit()
            return store

        def text(fraction):
            return '' if fraction is None else number_text(units.fraction_to_percent(fraction))

        first_field = self._field(f'carbon.{phase}', f'{name} carbon [%]', text(first), parser('first'),
                                  storer('first'), field_label='%' if choice != CARBON_RANGE else 'min %',
                                  classes='w-24')
        second_field = self._field(f'carbon.{phase}.max', f'{name} carbon maximum [%]', text(second),
                                   parser('second'), storer('second'), field_label='max %', classes='w-24')

        def show():
            first_field.set_visibility(select.value != 'free')
            second_field.set_visibility(select.value == CARBON_RANGE)

        def change(_):
            show()
            self.apply(f'carbon.{phase}.max')
            self.apply(f'carbon.{phase}')
        select.on_value_change(change)
        show()

    # ------------------------------------------------------------ aqueous phase
    def _aqueous(self):
        data = self.data
        ui.checkbox('Include an aqueous phase', value=get(data, ['aqueous', 'enabled'], True) is not False,
                    on_change=lambda e: self.set(['aqueous', 'enabled'], None if e.value else False)) \
            .mark('case-aqueous-enabled')
        ui.select(AQUEOUS_MODELS, label='Model', value=get(data, ['aqueous', 'model'], 'electrolyte'),
                  on_change=lambda e: self.set(['aqueous', 'model'], e.value)).classes('w-full')

        def store(value):
            put(self.data, ['aqueous', 'cutoff_C'], None if value is None else tidy(units.to_file(value, self.t_unit)))
        self.cutoff_field = self._field('aqueous.cutoff', units.field_label('Remove at and above', self.t_unit),
                                        format_numbers(get(data, ['aqueous', 'cutoff_C']), self.t_unit),
                                        lambda text: parse_number(text), store, classes='w-full')
        ui.label(CUTOFF_CAPTION).classes(CAPTION)
        self.conflict_label = ui.label().classes('text-body2 text-warning').mark('case-conflict')

    # ------------------------------------------------------------ solver
    def _solver(self):
        data = self.data
        value = get(data, ['solver', 'method'], 'ipopt')
        options = {**METHODS, **FD_METHOD} if value in FD_METHOD else METHODS
        method = ui.select(options, label='Method', value=value,
                           on_change=lambda e: self._set_method(e.value)).classes('w-full').mark('case-method')
        with method:
            ui.tooltip(METHOD_TOOLTIP).classes('text-body2 max-w-sm')
        starts = get(data, ['solver', 'starts'], 1)

        def store(value):
            put(self.data, ['solver', 'starts'], value)
            self.starts_tooltip.set_text(initial_guesses_tooltip(value))
        field = self._field('solver.starts', 'Initial guesses', number_text(starts),
                            lambda text: parse_number(text, required=True, minimum=1, integer=True), store,
                            classes='w-full')
        with field:
            self.starts_tooltip = ui.tooltip(initial_guesses_tooltip(starts)).classes('text-body2 max-w-sm') \
                .mark('case-starts-tooltip')
        continuation = ui.checkbox('Continuation', value=get(data, ['sweep', 'continuation'], True) is not False,
                                   on_change=lambda e: self.set(['sweep', 'continuation'], None if e.value else False))
        with continuation:
            ui.tooltip(CONTINUATION_TOOLTIP).classes('text-body2 max-w-sm')

    def _set_method(self, method):
        if getattr(self, 'fd_box', None) is not None:
            self.fd_box.set_visibility(method == 'slsqp-fd')
        self.set(['solver', 'method'], method)

    # ------------------------------------------------------------ advanced settings
    def _advanced(self):
        data = self.data
        with ui.expansion('Advanced settings', icon='tune', caption=ADVANCED_CAPTION).classes('w-full'):
            with ui.row().classes('w-full gap-8 items-start no-wrap'):
                with ui.column().classes('gap-1'):
                    ui.label('Starting point').classes('text-weight-medium')
                    ui.select(INITIAL_COMPOSITIONS, label='Initial composition',
                              value=get(data, ['initial_guess', 'initialization'], 'auto'),
                              on_change=lambda e: self.set(['initial_guess', 'initialization'],
                                                           None if e.value == 'auto' else e.value)).classes('w-72')
                    ui.label(INITIAL_COMPOSITION_CAPTION).classes(CAPTION + ' w-72')
                    self._optional('advanced.gas_water', 'Water in the gas at the start [fraction]',
                                   ['initial_guess', 'initial_gas_water_fraction'], '1e-6', minimum=0, maximum=1)
                    self._optional('advanced.water_floor', 'Smallest amount of aqueous water [mol / mol C]',
                                   ['aqueous', 'water_floor'], '1e-8', minimum=0)
                with ui.column().classes('gap-1'):
                    ui.label('Solver details').classes('text-weight-medium')
                    self._optional('advanced.max_iterations', 'Maximum iterations', ['solver', 'max_iterations'],
                                   '5000', integer=True, minimum=1)
                    self._optional('advanced.seed', 'Random seed for starting points', ['solver', 'seed'], '0',
                                   SEED_CAPTION, integer=True, minimum=0)
                    stability = {True: 'on', False: 'off'}.get(get(data, ['solver', 'stability_test']), 'default')
                    ui.select(STABILITY_TEST, label='Phase stability test', value=stability,
                              on_change=lambda e: self.set(['solver', 'stability_test'],
                                                           {'on': True, 'off': False}.get(e.value))).classes('w-72')
                    self._optional('advanced.energy_scale', 'Energy scale [J / mol C]', ['solver', 'energy_scale'],
                                   'automatic', ENERGY_SCALE_CAPTION, above=0)
                    self.fd_box = ui.column().classes('gap-1')
                    with self.fd_box:
                        self._optional('advanced.fd_ftol', 'Finite-difference tolerance', ['solver', 'fd_ftol'],
                                       '1e-9', above=0)
                        self._optional('advanced.fd_eps', 'Finite-difference step', ['solver', 'fd_eps'],
                                       '1e-7', above=0)
                    self.fd_box.set_visibility(get(data, ['solver', 'method'], 'ipopt') == 'slsqp-fd')
                with ui.column().classes('gap-1'):
                    ui.label('When are two results different?').classes('text-weight-medium')
                    ui.label(DISTINCT_CAPTION).classes(CAPTION + ' w-72')
                    self._optional('advanced.composition_tol', 'Composition [mol / mol C]',
                                   ['solver', 'distinct_composition_tol'], '1e-4', COMPOSITION_TOL_CAPTION,
                                   minimum=0)
                    self._optional('advanced.gibbs_tol', 'Gibbs energy [J / mol C]',
                                   ['solver', 'distinct_gibbs_tol'], '1', GIBBS_TOL_CAPTION, minimum=0)
            with ui.row().classes('w-full gap-8 items-start no-wrap q-mt-md'):
                with ui.column().classes('gap-1'):
                    ui.label('Minimum amounts [mol / mol C]').classes('text-weight-medium')
                    minimum = get(data, ['constraints', 'minimum_amounts']) or {}
                    self.minimum_items = [[p, s, v] for p, names in minimum.items() for s, v in names.items()]
                    self._minimum_rows = Part(self._build_minimum_rows)
                    ui.button('+ Add limit', on_click=self._add_limit).props('flat dense no-caps') \
                        .mark('case-add-limit')
                with ui.column().classes('gap-1'):
                    ui.label('Atomic weights').classes('text-weight-medium')
                    ui.label(ATOMIC_WEIGHTS_CAPTION).classes(CAPTION + ' w-72')
                    self.weight_items = [[e, v] for e, v in (data.get('atomic_weights') or {}).items()]
                    self._weight_rows = Part(self._build_weight_rows)
                    with ui.row().classes('items-start gap-2'):
                        symbol = ui.input('Element symbol').classes('w-32').mark('case-new-weight')
                        ui.button('+ Add element', on_click=lambda: self._add_weight(symbol)) \
                            .props('flat dense no-caps').classes('q-mt-sm').mark('case-add-weight')
                with ui.column().classes('gap-1'):
                    ui.label('Measured data to show in plots').classes('text-weight-medium')
                    self._measured()

    # minimum amounts
    def _build_minimum_rows(self):
        self._forget('advanced.min')
        if self.minimum_items:
            with ui.grid(columns='9rem 14rem 7rem 2rem').classes('items-start gap-x-2 gap-y-0'):
                for text in ('Phase', 'Species', 'Minimum', ''):
                    ui.label(text).classes(CAPTION)
                for i, (phase, species, value) in enumerate(self.minimum_items):
                    ui.select(PHASE_NAMES, value=phase if phase in PHASE_NAMES else None,
                              on_change=lambda e, i=i: self._limit_part(i, 0, e.value)).classes('w-36')
                    options = self.state.phase_options.get(phase, []) + ([species] if species else [])
                    select = ui.select(list(dict.fromkeys(options)), value=species or None, with_input=True,
                                       on_change=lambda e, i=i: self._limit_part(i, 1, e.value)).classes('w-56')
                    self.fields[f'advanced.min-species:{i}'] = (select.without_auto_validation(), 'Minimum amounts',
                                                                 None, None)
                    self._field(f'advanced.min:{i}', 'Minimum amounts', number_text(value),
                                lambda text: parse_number(text, required=True, minimum=0),
                                lambda v, i=i: self._limit_part(i, 2, v, refresh=False), field_label='',
                                classes='w-28')
                    ui.button(icon='close', on_click=lambda i=i: self._remove_limit(i)) \
                        .props('flat dense round size=sm').classes('q-mt-sm')
        self._commit_limits()

    def _limit_part(self, i, part, value, refresh=True):
        self.minimum_items[i][part] = value
        if part == 0:
            self.minimum_items[i][1] = ''
            self._minimum_rows.refresh()
        self._commit_limits()
        if refresh:
            self.changed()

    def _add_limit(self):
        self.minimum_items.append(['gas', '', 1e-8])
        self._minimum_rows.refresh()
        self.changed()

    def _remove_limit(self, i):
        self.minimum_items.pop(i)
        self._minimum_rows.refresh()
        self.changed()

    def _commit_limits(self):
        out = {}
        for i, (phase, species, value) in enumerate(self.minimum_items):
            key = f'advanced.min-species:{i}'
            if not phase or not species:
                self.parse_errors[key] = 'Choose a phase and a species, or remove the row.'
            elif species in out.get(phase, {}):
                self.parse_errors[key] = f'{species} ({PHASE_NAMES[phase].lower()}) is listed twice.'
            else:
                self.parse_errors.pop(key, None)
                if value is not None:
                    out.setdefault(phase, {})[species] = value
        put(self.data, ['constraints', 'minimum_amounts'], out or None)

    # atomic weights
    def _build_weight_rows(self):
        self._forget('advanced.aw')
        if self.weight_items:
            with ui.grid(columns='3rem 8rem 2rem').classes('items-start gap-x-2 gap-y-0'):
                for text in ('Element', 'g/mol', ''):
                    ui.label(text).classes(CAPTION)
                for i, (element, value) in enumerate(self.weight_items):
                    ui.label(element).classes('q-mt-sm')
                    self._field(f'advanced.aw:{i}', f'Atomic weight of {element}', number_text(value),
                                lambda text: parse_number(text, required=True, above=0),
                                lambda v, i=i: self._weight_value(i, v), field_label='', classes='w-28')
                    ui.button(icon='close', on_click=lambda i=i: self._remove_weight(i)) \
                        .props('flat dense round size=sm').classes('q-mt-sm')

    def _commit_weights(self):
        weights = {e: v for e, v in self.weight_items if v is not None}
        put(self.data, ['atomic_weights'], weights or None)

    def _weight_value(self, i, value):
        self.weight_items[i][1] = value
        self._commit_weights()

    def _add_weight(self, symbol_input):
        symbol = element_symbol(symbol_input.value)
        error = element_error(symbol, [e for e, _ in self.weight_items])
        symbol_input.error = error
        if error:
            return
        self.weight_items.append([symbol, ATOMIC_WEIGHTS.get(symbol)])
        self._commit_weights()
        self._weight_rows.refresh()
        self.changed()

    def _remove_weight(self, i):
        self.weight_items.pop(i)
        self._commit_weights()
        self._weight_rows.refresh()
        self.changed()

    # measured data
    def _measured(self):
        data = self.data
        files = self.state.comparison_files(self.document.base_dir())
        current = get(data, ['comparison', 'file'])
        options = {'': '(none)', **files, **({current: current} if current and current not in files else {})}
        file = ui.select(options, label='CSV file', value=current or '').classes('w-80')
        quantity = ui.select(QUANTITY_CHOICES, label='Its unit',
                             value=get(data, ['comparison', 'quantity'], 'mol_per_kg_dry_feed')).classes('w-80')
        label = ui.input('Label in the legend', value=get(data, ['comparison', 'label'], '')).classes('w-80')
        ui.label(MEASURED_CAPTION).classes(CAPTION + ' w-80')

        def commit(_):
            if not file.value:
                comparison = None
            else:
                comparison = {'file': file.value, 'quantity': quantity.value}
                if (label.value or '').strip():
                    comparison['label'] = label.value
            self.set(['comparison'], comparison)
        for element in (file, quantity, label):
            element.on_value_change(commit)

    # ------------------------------------------------------------ derived numbers
    def _update_derived(self):
        data = self.data
        self.point_label.set_text(point_count_text(point_count(data.get('conditions'),
                                                               get(data, ['aqueous', 'cutoff_C']))))
        self._update_feed()
        notice = conflict_notice(data, self.t_unit)
        self.conflict_label.set_text(notice or '')
        self.conflict_label.set_visibility(bool(notice))
        for key, (element, _, _, _) in self.fields.items():
            element.error = self.parse_errors.get(key) or self.feed_errors.get(key)


# ---------------------------------------------------------------- small UI helpers
class Part:
    """A part of the page that ``refresh()`` rebuilds at once.

    (``ui.refreshable`` rebuilds in a background task; the form needs the new
    elements right away to fill in the derived numbers.)
    """

    def __init__(self, build, classes='w-full gap-2'):
        self.build = build
        self.box = ui.column().classes(classes)
        self.refresh()

    def refresh(self):
        self.box.clear()
        with self.box:
            self.build()


def card(title, caption=None, grow=False):
    """A titled card (use as a context manager)."""
    element = ui.card().classes('grow basis-0' if grow else 'w-full')
    with element:
        ui.label(title).classes(TITLE)
        if caption:
            ui.label(caption).classes(CAPTION)
    return element


def set_hint(element, text):
    element._props['hint'] = text
    element.update()
