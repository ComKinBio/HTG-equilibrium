"""The new / edit species dialog of the Species tab: a form over a SpeciesSpec and the live check.

    SpeciesDialog(model, user_dir, on_saved).open()                     # new species
    SpeciesDialog(model, user_dir, on_saved, spec=spec_from_model(model, 'methane'),
                  editing='methane', overridden=False).open()            # edit a species

on_saved(name, action) is called after the user database changed (action 'saved'
or 'reset'); it may be a coroutine function.
"""
import copy
import inspect
import math

from nicegui import ui

from ..formulas import formula_summary, hill_formula
from ..saft import normalize_ion_name
from ..user_species import (ASSOCIATION_SCHEMES, FIELD_LABELS, HKF_FIELDS, KIND_LABELS, KINDS,
                            MAIER_KELLEY_FIELDS, SUPCRT_UNITS, SpeciesSpec, UserDatabase, check_species,
                            find_supcrt_node, nasa_match_message, nasa_records_for_formula,
                            parse_nasa_record, supcrt_aqueous_records)

THERMO_OPTIONS = {
    'molecule': {'nasa_record': 'Use a NASA-7 record already in the database',
                 'nasa': 'Enter NASA-7 coefficients',
                 'paste': 'Paste a NASA-7 record (4 lines, CHEMKIN format)'},
    'ion': {'hkf_record': 'Use an aqueous HKF record from SUPCRT92',
            'hkf': 'Enter HKF parameters (ΔGf°, ΔHf°, S°, a1–a4, c1, c2, ω)'},
    'solid': {'nasa': 'NASA-7 coefficients (phase S)',
              'maier_kelley': 'Maier–Kelley: ΔGf°, ΔHf°, S°, V°, a, b, c and maximum temperature',
              'nasa_record': 'Use a NASA-7 record already in the database'},
}
SOLID_NOTE = 'Single temperature range only. Minerals with phase transitions are not supported yet.'
SAVE_NOTE = 'Saved to your own database. The shipped databases and example cases are not changed.'
EOS_MODELS = ['PC-SAFT', 'ePC-SAFT']
NOT_A_NUMBER = 'Enter a number'


def parse_number(text):
    """The number in a text field; None when empty or not a number (1.5E-3 and 1.5D-3 both work)."""
    text = str(text if text is not None else '').strip().replace('D', 'E').replace('d', 'e')
    try:
        value = float(text)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def format_number(value):
    """Text for a number field: short where that is exact, else full precision; '' for none."""
    try:
        value = float(value)
    except (TypeError, ValueError):
        return ''
    if not math.isfinite(value):
        return ''
    short = f'{value:.12g}'
    return short if float(short) == value else repr(value)


def _number_ok(text):
    return None if not str(text or '').strip() or parse_number(text) is not None else NOT_A_NUMBER


def _field_label(key):
    label = FIELD_LABELS.get(key, key)
    return f'{label[:1].upper()}{label[1:]} [{SUPCRT_UNITS[key]}]'


def _same_name(a, b):
    return normalize_ion_name(str(a or '').strip()) == normalize_ion_name(str(b or '').strip())


class SpeciesDialog:
    """The dialog (see the module docstring).

    model       the model in use (records to choose from, the name check)
    user_dir    the user database folder written on save
    editing     name of the species being edited (None: new species)
    overridden  the edited species is a shipped one with user values: show "Reset to original"
    user_only   the edited species exists only in the user database: renaming it moves it
    """

    def __init__(self, model, user_dir, on_saved, spec=None, editing=None, overridden=False,
                 user_only=False):
        self.model = model
        self.user_dir = user_dir
        self.on_saved = on_saved
        self.editing = editing
        self.user_only = user_only
        s = self.spec = copy.deepcopy(spec) if spec is not None else SpeciesSpec()
        nasa = dict(s.nasa or {})
        s.nasa = {key: nasa.get(key) for key in ('T_low', 'T_mid', 'T_high')}
        for key in ('high', 'low'):
            values = list(nasa.get(key) or [])
            s.nasa[key] = (values + [None] * 7)[:7]
        s.hkf = {**{k: None for k in HKF_FIELDS}, **(s.hkf or {})}
        s.maier_kelley = {**{k: None for k in MAIER_KELLEY_FIELDS}, **(s.maier_kelley or {})}
        self.kij = [[other, kij] for other, kij in s.kij]
        options = THERMO_OPTIONS[s.kind]
        self.mode = s.thermo if s.thermo in options else next(iter(options))
        s.thermo = 'nasa' if self.mode == 'paste' else self.mode
        self.inputs = {}
        self.molecules = sorted(model.registry.all_molecular, key=str.casefold)

        self.dialog = ui.dialog().props('persistent')
        with self.dialog, ui.card().classes('w-[80rem] max-w-full'):
            ui.label(f'Edit {editing}' if editing else 'New species').classes('text-h6')
            with ui.row().classes('w-full no-wrap items-start gap-6'):
                with ui.column().classes('grow min-w-0 gap-4 max-h-[70vh] overflow-auto pr-2'):
                    self._kind()
                    self._name()
                    self._thermo()
                    self._eos()
                    self._source()
                with ui.column().classes('w-80 shrink-0 gap-2'):
                    ui.label('Check').classes('text-subtitle1 text-weight-medium')
                    self.check_box = ui.column().classes('gap-1')
                    self.can_be_in = ui.label().classes('text-body2')
                    ui.label(SAVE_NOTE).classes('text-caption text-grey-7')
            with ui.row().classes('w-full items-center'):
                if overridden:
                    ui.button('Reset to original', on_click=self.reset).props('flat no-caps color=negative')
                ui.space()
                ui.button('Cancel', on_click=self.dialog.close).props('flat no-caps')
                ui.button('Save species', on_click=self.save).props('no-caps')
        self.refresh_records()
        self.show_parts()
        self.update_check()

    def open(self):
        self.dialog.open()
        return self

    # ------------------------------------------------------------ helpers
    def _number(self, label, value, on_change, classes='w-40', mark=None):
        field = ui.input(label, value=format_number(value), validation=_number_ok,
                         on_change=lambda e: (on_change(parse_number(e.value)), self.update_check()))
        if mark:
            field.mark(mark)
        return field.props('dense').classes(classes)

    def _set(self, key, value):
        setattr(self.spec, key, value)
        self.update_check()

    @staticmethod
    def _section(title):
        card = ui.column().classes('w-full gap-2')
        with card:
            ui.label(title).classes('text-subtitle1 text-weight-medium')
        return card

    # ------------------------------------------------------------ sections
    def _kind(self):
        with ui.row().classes('items-center gap-3'):
            ui.label('Kind').classes('text-subtitle1 text-weight-medium')
            self.kind_toggle = ui.toggle({k: KIND_LABELS[k] for k in KINDS}, value=self.spec.kind,
                                         on_change=self.set_kind).props('dense no-caps')
            if self.editing:
                self.kind_toggle.disable()

    def _name(self):
        s = self.spec
        with self._section('Name and formula'):
            with ui.row().classes('items-start gap-3'):
                ui.input('Name', value=s.name, on_change=lambda e: self._set('name', e.value)) \
                    .props('dense').classes('w-56').mark('species-name')
                self.formula_input = ui.input('Formula', value=s.formula, on_change=self.set_formula) \
                    .props('dense').classes('w-48').mark('species-formula')
                self.charge = ui.number('Charge', value=s.charge, step=1, precision=0,
                                        on_change=self.set_charge).props('dense').classes('w-24')
            self.summary = ui.label().classes('text-body2 text-grey-8')

    def _thermo(self):
        s = self.spec
        with self._section('Thermochemistry'):
            self.thermo_radio = ui.radio(THERMO_OPTIONS[s.kind], value=self.mode, on_change=self.set_mode) \
                .props('dense')
            self.nasa_record_box = ui.column().classes('gap-1')
            with self.nasa_record_box:
                self.nasa_select = ui.select([], label='NASA-7 record', with_input=True, clearable=True,
                                             on_change=self.set_nasa_record).props('dense').classes('w-64')
                self.nasa_match = ui.label().classes('text-caption text-grey-8')
            self.nasa_box = ui.column().classes('gap-1')
            with self.nasa_box:
                with ui.row().classes('gap-3'):
                    for key, label in (('T_low', 'T low [K]'), ('T_mid', 'T common [K]'),
                                       ('T_high', 'T high [K]')):
                        self.inputs[key] = self._number(label, s.nasa[key],
                                                        lambda v, k=key: s.nasa.__setitem__(k, v), 'w-32')
                for key, title in (('high', 'High-range coefficients (T common to T high)'),
                                   ('low', 'Low-range coefficients (T low to T common)')):
                    ui.label(title).classes('text-body2')
                    with ui.row().classes('gap-2'):
                        for i in range(7):
                            self.inputs[f'{key}{i}'] = self._number(
                                f'a{i + 1}', s.nasa[key][i],
                                lambda v, k=key, j=i: s.nasa[k].__setitem__(j, v), 'w-28')
            self.paste_box = ui.column().classes('w-full gap-1')
            with self.paste_box:
                ui.textarea('NASA-7 record', on_change=self.set_paste) \
                    .props('dense outlined input-class="font-mono text-xs" rows=5').classes('w-full') \
                    .mark('nasa-paste')
                self.paste_message = ui.label().classes('text-caption text-grey-8')
            self.hkf_record_box = ui.column().classes('gap-1')
            with self.hkf_record_box:
                self.hkf_select = ui.select([], label='SUPCRT92 record', with_input=True, clearable=True,
                                            on_change=self.set_hkf_record).props('dense').classes('w-64')
                self.hkf_match = ui.label().classes('text-caption text-grey-8')
            self.hkf_box = ui.row().classes('gap-2')
            with self.hkf_box:
                for key in HKF_FIELDS:
                    self._number(_field_label(key), s.hkf.get(key),
                                 lambda v, k=key: s.hkf.__setitem__(k, v), 'w-44')
            self.mk_box = ui.column().classes('gap-1')
            with self.mk_box:
                with ui.row().classes('gap-2'):
                    for key in MAIER_KELLEY_FIELDS:
                        self._number(_field_label(key), s.maier_kelley.get(key),
                                     lambda v, k=key: s.maier_kelley.__setitem__(k, v), 'w-44')
                ui.label(SOLID_NOTE).classes('text-caption text-grey-8')

    def _eos(self):
        s = self.spec
        self.eos_box = self._section('Equation of state')
        with self.eos_box:
            with ui.row().classes('items-start gap-3'):
                self.model_select = ui.select(EOS_MODELS, value=s.model if s.model in EOS_MODELS else 'PC-SAFT',
                                              label='Model', on_change=lambda e: self._set('model', e.value)) \
                    .props('dense').classes('w-32')
                self._number('Segment number m', s.m, lambda v: setattr(s, 'm', v), mark='eos-m')
                self._number('Segment diameter σ [Å]', s.sigma, lambda v: setattr(s, 'sigma', v), mark='eos-sigma')
                self._number('Energy ε/k [K]', s.epsilon, lambda v: setattr(s, 'epsilon', v),
                             mark='eos-epsilon')
            self.assoc_box = ui.row().classes('items-start gap-3')
            with self.assoc_box:
                scheme = s.association_scheme if s.association_scheme in ASSOCIATION_SCHEMES else '-'
                ui.select(ASSOCIATION_SCHEMES, value=scheme, label='Association scheme',
                          on_change=lambda e: self._set('association_scheme', e.value)) \
                    .props('dense').classes('w-40')
                self._number('Association energy εAB/k [K]', s.epsilonAB, lambda v: setattr(s, 'epsilonAB', v),
                             'w-52')
                self._number('Association volume κAB', s.kappaAB, lambda v: setattr(s, 'kappaAB', v), 'w-44')
            self.kij_box = ui.column().classes('gap-1')
            with self.kij_box:
                self.kij_rows()
                ui.button('+ Binary interaction (kij)', on_click=self.add_kij).props('flat dense no-caps')

    @ui.refreshable_method
    def kij_rows(self):
        for row in self.kij:
            with ui.row().classes('items-center gap-3'):
                options = list(dict.fromkeys(self.molecules + ([row[0]] if row[0] else [])))
                ui.select(options, value=row[0] or None, label='Other species', with_input=True,
                          on_change=lambda e, r=row: (r.__setitem__(0, e.value or ''), self.update_check())) \
                    .props('dense').classes('w-56')
                ui.input('kij', value=format_number(row[1]), validation=_number_ok,
                         on_change=lambda e, r=row: (r.__setitem__(1, parse_number(e.value)),
                                                     self.update_check())).props('dense').classes('w-28')
                ui.button(icon='close', on_click=lambda r=row: self.remove_kij(r)) \
                    .props('flat dense round size=sm').tooltip('Remove this kij')

    def _source(self):
        s = self.spec
        with self._section('Source'):
            with ui.row().classes('w-full items-start gap-3'):
                ui.input('Authors', value=s.authors, on_change=lambda e: self._set('authors', e.value)) \
                    .props('dense').classes('grow')
                ui.input('Year', value=s.year, on_change=lambda e: self._set('year', e.value)) \
                    .props('dense').classes('w-24')
                ui.input('DOI', value=s.doi, on_change=lambda e: self._set('doi', e.value)) \
                    .props('dense').classes('w-64')

    # ------------------------------------------------------------ events
    def set_kind(self, e):
        s = self.spec
        s.kind = e.value
        options = THERMO_OPTIONS[s.kind]
        if self.mode not in options:
            self.mode = next(iter(options))
        s.thermo = 'nasa' if self.mode == 'paste' else self.mode
        self.thermo_radio.set_options(options, value=self.mode)
        if s.kind == 'ion':
            s.model = 'ePC-SAFT'
            self.model_select.set_value('ePC-SAFT')
        self.refresh_records()
        self.show_parts()
        self.update_check()

    def set_mode(self, e):
        if not e.value:
            return
        self.mode = e.value
        self.spec.thermo = 'nasa' if self.mode == 'paste' else self.mode
        self.show_parts()
        self.update_check()

    def set_formula(self, e):
        self.spec.formula = e.value or ''
        self.refresh_records()
        self.update_check()

    def set_charge(self, e):
        try:
            self.spec.charge = int(e.value or 0)
        except (TypeError, ValueError):
            self.spec.charge = 0
        self.refresh_records()
        self.update_check()

    def set_nasa_record(self, e):
        self.spec.nasa_record = e.value or ''
        self._nasa_match()
        self.update_check()

    def set_hkf_record(self, e):
        self.spec.supcrt_record = e.value or ''
        self._hkf_match()
        self.update_check()

    def set_paste(self, e):
        text = e.value or ''
        if not text.strip():
            self.paste_message.set_text('')
            return
        try:
            record = parse_nasa_record(text)
        except ValueError as exc:
            self.paste_message.set_text(str(exc))
            self.update_check()
            return
        s = self.spec
        for key in ('T_low', 'T_mid', 'T_high'):
            self.inputs[key].set_value(format_number(record[key]))
        for key in ('high', 'low'):
            for i, value in enumerate(record[key]):
                self.inputs[f'{key}{i}'].set_value(format_number(value))
        formula = hill_formula(record['elements']) if record['elements'] else ''
        if not s.formula and formula:
            self.formula_input.set_value(formula)
        self.paste_message.set_text(f'Read record {record["name"]}: {formula}, '
                                    f'{record["T_low"]:g}–{record["T_high"]:g} K')
        self.update_check()

    def add_kij(self):
        self.kij.append(['', None])
        self.kij_rows.refresh()
        self.update_check()

    def remove_kij(self, row):
        self.kij.remove(row)
        self.kij_rows.refresh()
        self.update_check()

    # ------------------------------------------------------------ display
    def show_parts(self):
        kind, mode = self.spec.kind, self.mode
        self.charge.set_visibility(kind == 'ion')
        self.nasa_record_box.set_visibility(mode == 'nasa_record')
        self.nasa_box.set_visibility(mode == 'nasa')
        self.paste_box.set_visibility(mode == 'paste')
        self.hkf_record_box.set_visibility(mode == 'hkf_record')
        self.hkf_box.set_visibility(mode == 'hkf')
        self.mk_box.set_visibility(mode == 'maier_kelley')
        self.eos_box.set_visibility(kind != 'solid')
        self.model_select.set_enabled(kind != 'ion')
        self.assoc_box.set_visibility(kind == 'molecule')
        self.kij_box.set_visibility(kind == 'molecule')

    def refresh_records(self):
        """Record choices that follow the formula, kind and charge."""
        s = self.spec
        solid = s.kind == 'solid'
        names = nasa_records_for_formula(self.model, s.formula, solid=solid) if s.formula else []
        if not s.nasa_record and len(names) == 1 and self.mode == 'nasa_record':
            s.nasa_record = names[0]
        options = list(dict.fromkeys(names + ([s.nasa_record] if s.nasa_record else [])))
        self.nasa_select.set_options(options, value=s.nasa_record or None)
        self._nasa_match()
        if s.kind == 'ion':
            names = supcrt_aqueous_records(self.model, s.charge if s.charge else None)
            options = list(dict.fromkeys(names + ([s.supcrt_record] if s.supcrt_record else [])))
            self.hkf_select.set_options(options, value=s.supcrt_record or None)
            self._hkf_match()

    def _nasa_match(self):
        s = self.spec
        if not s.formula:
            self.nasa_match.set_text('Enter the formula to list the records.')
            return
        self.nasa_match.set_text(nasa_match_message(self.model, s.formula, s.nasa_record or None,
                                                    solid=s.kind == 'solid'))

    def _hkf_match(self):
        s = self.spec
        node = find_supcrt_node(self.model.databases.supcrt, s.supcrt_record) if s.supcrt_record else None
        if node is None:
            self.hkf_match.set_text('')
            return
        charge = int(float(node.findtext('Charge') or 0))
        self.hkf_match.set_text('✓ Charge matches' if charge == s.charge else
                                f'! Charge of {s.supcrt_record} is {charge:+d}, not {s.charge:+d}')

    def current_spec(self):
        spec = copy.deepcopy(self.spec)
        spec.thermo = 'nasa' if self.mode == 'paste' else self.mode
        spec.kij = [(other, kij) for other, kij in self.kij if other or kij is not None]
        return spec

    def check(self):
        return check_species(self.current_spec(), self.model, editing=self.editing)

    def update_check(self):
        s = self.spec
        try:
            if s.formula:
                self.summary.set_text(formula_summary(s.formula, s.charge if s.kind == 'ion' else None))
            else:
                self.summary.set_text('')
        except ValueError as exc:
            self.summary.set_text(str(exc))
        result = self.check()
        self.check_box.clear()
        with self.check_box:
            for item in result.items:
                ui.label(f'{"✓" if item.ok else "!"} {item.text}') \
                    .classes('text-body2 ' + ('text-positive' if item.ok else 'text-orange-9'))
        self.can_be_in.set_text(result.can_be_in_text)

    # ------------------------------------------------------------ save
    async def save(self):
        spec = self.current_spec()
        try:
            db = UserDatabase.load(self.user_dir)
            name = db.save_species(spec, self.model)
            if self.user_only and self.editing and not _same_name(name, self.editing):
                db.remove_species(self.editing)
        except (ValueError, OSError) as exc:
            ui.notify(f'Not saved: {exc}', type='negative')
            return
        self.dialog.close()
        await _call(self.on_saved, name, 'saved')

    async def reset(self):
        try:
            UserDatabase.load(self.user_dir).reset_species(self.editing)
        except OSError as exc:
            ui.notify(f'Not reset: {exc}', type='negative')
            return
        self.dialog.close()
        await _call(self.on_saved, self.editing, 'reset')


async def _call(function, *args):
    result = function(*args)
    if inspect.isawaitable(result):
        await result

