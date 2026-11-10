"""The case editor: a form over a case file (CaseDocument.data)."""
import re
from pathlib import Path

from nicegui import run as background
from nicegui import ui

from ..case import INITIALIZATIONS, SOLVERS
from .document import CaseDocument

PHASE_NAMES = {'gas': 'Gas', 'aqueous': 'Aqueous', 'organic_liquid': 'Organic liquid', 'solid': 'Solid'}
QUANTITY_CHOICES = {'mol_per_kg_dry_feed': 'mol / kg dry feed', 'mol_per_mol_C': 'mol / mol dry-feed C',
                    'g_per_kg_dry_feed': 'g / kg dry feed'}


# ---------------------------------------------------------------- data helpers
def get(data, path, default=None):
    for key in path:
        if not isinstance(data, dict) or key not in data:
            return default
        data = data[key]
    return data


def put(data, path, value):
    """Set data[path] = value; None (or '') removes the key."""
    *parents, last = path
    for key in parents:
        data = data.setdefault(key, {})
    if value is None or value == '':
        data.pop(last, None)
    else:
        data[last] = value


def parse_numbers(text):
    """'100, 150 200' or 'start:stop:step' (inclusive) -> list of floats."""
    text = (text or '').strip()
    if not text:
        return []
    values = []
    for part in re.split(r'[,;\s]+', text):
        if ':' in part:
            start, stop, step = (float(x) for x in part.split(':'))
            if step <= 0:
                raise ValueError('step must be positive')
            n = int(round((stop - start) / step))
            values += [round(start + i * step, 8) for i in range(n + 1) if start + i * step <= stop + 1e-9]
        elif part:
            values.append(float(part))
    return values


def format_numbers(value):
    if value is None:
        return ''
    values = value if isinstance(value, list) else [value]
    return ', '.join(f'{v:g}' for v in values)


def parse_points(text):
    """Lines 'T, P' -> [[T, P], ...]."""
    points = []
    for line in (text or '').splitlines():
        if line.strip():
            numbers = parse_numbers(line)
            if len(numbers) != 2:
                raise ValueError(f'"{line.strip()}": give temperature and pressure')
            points.append(numbers)
    return points


def parse_elements(text):
    """'Na=1, H=1 C=1' -> {'Na': 1, 'H': 1, 'C': 1}."""
    out = {}
    for part in re.split(r'[,;\s]+', (text or '').strip()):
        if part:
            if '=' not in part:
                raise ValueError(f'"{part}": write element=count, e.g. Na=1')
            element, count = part.split('=', 1)
            number = float(count)
            out[element.strip()] = int(number) if number.is_integer() else number
    return out


def format_elements(elements):
    return ', '.join(f'{e}={c:g}' for e, c in (elements or {}).items())


def scalar(value):
    return value[0] if isinstance(value, list) and len(value) == 1 else value


# ---------------------------------------------------------------- editor
class CaseEditor:
    """File controls, the form, a status line and the actions (check, queue)."""

    def __init__(self, state, on_queued):
        self.state = state
        self.on_queued = on_queued
        with ui.row().classes('w-full items-center gap-2'):
            self.file_select = ui.select(self._case_files(), label='Case file', with_input=True,
                                         on_change=self._open_selected).classes('w-80')
            ui.button('New', icon='note_add', on_click=self._new).props('flat')
            ui.button('Save', icon='save', on_click=self._save).props('flat')
            ui.button('Save as', icon='save_as', on_click=self._save_as).props('flat')
            ui.space()
            ui.button('Check', icon='fact_check', on_click=self._check).props('outline')
            ui.button('Add to queue', icon='play_arrow', on_click=self._queue).props('color=primary')
        self.file_label = ui.label().classes('text-caption text-grey-7')
        self.status = ui.label().classes('text-body2')
        self.form()
        self.changed()

    # ------------------------------------------------------------ files
    @property
    def document(self):
        return self.state.document

    @property
    def data(self):
        return self.state.document.data

    def _case_files(self):
        root = self.state.cases_dir
        files = sorted(root.glob('**/*.toml')) if root.exists() else []
        return {str(p): str(p.relative_to(root)) for p in files}

    def _open_selected(self, event):
        if event.value:
            self.open(event.value)

    def open(self, path):
        try:
            self.state.document = CaseDocument.open(path)
        except Exception as exc:
            ui.notify(f'Could not open {path}: {exc}', type='negative', multi_line=True)
            return
        self.form.refresh()
        self.changed()

    def _new(self):
        self.state.document = CaseDocument.new()
        self.file_select.set_value(None)
        self.form.refresh()
        self.changed()

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
        try:
            saved = self.document.save(path)
        except Exception as exc:
            ui.notify(f'Not saved: {exc}', type='negative', multi_line=True)
            return
        ui.notify(f'Saved {saved}', type='positive')
        self.file_select.set_options(self._case_files())
        self.changed()

    # ------------------------------------------------------------ status and actions
    def set(self, path, value):
        put(self.data, path, value)
        self.changed()

    def changed(self):
        where = str(self.document.path) if self.document.path else 'not saved yet'
        self.file_label.set_text(f'File: {where}')
        try:
            case = self.document.case()
            points, pieces = len(case.grid()), len(case.chains())
            self.status.set_text(f'✓ {points} point(s) in {pieces} independent piece(s)'
                                 + ('' if pieces > 1 or points == 1 else
                                    ' — continuation runs them in order on one worker'))
            self.status.classes(replace='text-body2 text-positive')
        except Exception as exc:
            self.status.set_text(f'✗ {exc}')
            self.status.classes(replace='text-body2 text-negative')

    async def _check(self):
        try:
            case = self.document.case()
        except Exception as exc:
            ui.notify(str(exc), type='negative', multi_line=True)
            return
        notification = ui.notification('Checking every point…', spinner=True, timeout=None)
        try:
            table = await background.io_bound(self.state.model.check, case)
        finally:
            notification.dismiss()
        bad = table[table.status != 'ready']
        with ui.dialog() as dialog, ui.card().classes('w-[56rem] max-w-full'):
            if bad.empty:
                ui.label(f'All {len(table)} points passed the checks.').classes('text-positive text-h6')
            else:
                ui.label(f'{len(bad)} of {len(table)} points will probably fail:').classes('text-negative text-h6')
                rows = [{'row': i, **r} for i, r in enumerate(bad[['T_C', 'P_bar', 'error']].to_dict('records'))]
                ui.table(rows=rows, columns=[{'name': c, 'label': c, 'field': c, 'align': 'left'}
                                             for c in ('T_C', 'P_bar', 'error')],
                         row_key='row', pagination=15).classes('w-full')
            ui.label('The checks need no equation of state; passing does not guarantee convergence.') \
                .classes('text-caption')
            ui.button('Close', on_click=dialog.close).props('flat')
        dialog.open()

    def _queue(self):
        try:
            case = self.document.case()
            run = self.state.submit(case)
        except Exception as exc:
            ui.notify(f'Not queued: {exc}', type='negative', multi_line=True)
            return
        ui.notify(f'Queued {case.name}: {run.total} point(s)', type='positive')
        self.on_queued(run)

    # ------------------------------------------------------------ the form
    @ui.refreshable_method
    def form(self):
        data = self.data
        with ui.column().classes('w-full gap-3'):
            with ui.card().classes('w-full'):
                with ui.row().classes('w-full'):
                    ui.input('Name (also the results folder)', value=data.get('name', ''),
                             on_change=lambda e: self.set(['name'], e.value)).classes('w-64')
                    ui.input('Description', value=data.get('description', ''),
                             on_change=lambda e: self.set(['description'], e.value)).classes('grow')
            self._conditions()
            self._feed()
            self._phases()
            with ui.row().classes('w-full items-stretch no-wrap gap-3'):
                self._aqueous()
                self._constraints()
                self._solver()
            self._advanced()

    def _section(self, title, caption=None):
        card = ui.card().classes('w-full')
        with card:
            ui.label(title).classes('text-subtitle1 text-weight-medium')
            if caption:
                ui.label(caption).classes('text-caption text-grey-7')
        return card

    # conditions
    def _conditions(self):
        data = self.data
        with self._section('Conditions', 'Temperatures in °C, pressures in bar (absolute). '
                                         'Lists: "300, 350, 400" or "300:450:25" (start:stop:step).'):
            explicit = get(data, ['conditions', 'points']) is not None
            mode = ui.toggle({'grid': 'Grid (every T at every P)', 'points': 'Explicit points'},
                             value='points' if explicit else 'grid').props('dense')
            grid_row = ui.row().classes('w-full')
            with grid_row:
                ui.input('Temperatures [°C]', value=format_numbers(get(data, ['conditions', 'temperature_C'])),
                         validation=self._numbers_validator(['conditions', 'temperature_C'])).classes('grow')
                ui.input('Pressures [bar]', value=format_numbers(get(data, ['conditions', 'pressure_bar'])),
                         validation=self._numbers_validator(['conditions', 'pressure_bar'])).classes('w-80')
            points_row = ui.row().classes('w-full')
            with points_row:
                text = '\n'.join(f'{t:g}, {p:g}' for t, p in (get(data, ['conditions', 'points']) or []))
                ui.textarea('Points: one "temperature, pressure" per line, solved in this order', value=text,
                            validation=self._points_validator).classes('w-full').props('autogrow')
            grid_row.bind_visibility_from(mode, 'value', value='grid')
            points_row.bind_visibility_from(mode, 'value', value='points')
            mode.on_value_change(lambda e: self._switch_conditions(e.value))

    def _switch_conditions(self, mode):
        conditions = self.data.setdefault('conditions', {})
        if mode == 'points' and 'points' not in conditions:
            try:
                grid = self.document.case().grid()
            except Exception:
                grid = []
            conditions.clear()
            conditions['points'] = [[t, p] for t, p in grid]
        elif mode == 'grid' and 'points' in conditions:
            points = conditions.pop('points') or [[300.0, 240.0]]
            conditions['temperature_C'] = sorted({t for t, _ in points})
            conditions['pressure_bar'] = scalar(sorted({p for _, p in points}))
        self.form.refresh()
        self.changed()

    def _numbers_validator(self, path):
        def validate(text):
            try:
                values = parse_numbers(text)
            except ValueError as exc:
                return str(exc)
            if not values:
                return 'Give at least one value'
            self.set(path, scalar(values))
            return None
        return validate

    def _points_validator(self, text):
        try:
            points = parse_points(text)
        except ValueError as exc:
            return str(exc)
        if not points:
            return 'Give at least one point'
        self.set(['conditions', 'points'], points)
        return None

    # feed
    def _feed(self):
        data = self.data
        feed = data.setdefault('feed', {})
        with self._section('Feed', 'Dry feed as an elemental analysis (mass %, dry basis) or as moles of '
                                   'species; water and salts per gram of dry feed.'):
            by_species = feed.get('species_moles') is not None
            mode = ui.toggle({'elements': 'Elemental analysis', 'species': 'Species moles'},
                             value='species' if by_species else 'elements').props('dense')
            with ui.row().classes('w-full items-start gap-8'):
                elements_box = ui.column()
                with elements_box:
                    self._mapping(['feed', 'elements_wt_percent'], 'Element', 'wt %')
                species_box = ui.column()
                with species_box:
                    self._mapping(['feed', 'species_moles'], 'Species', 'relative mol',
                                  key_options=self.state.feed_species)
                elements_box.bind_visibility_from(mode, 'value', value='elements')
                species_box.bind_visibility_from(mode, 'value', value='species')
                with ui.column():
                    ui.number('Added water [g / g dry feed]', value=feed.get('water_g_per_g_dry'), min=0,
                              on_change=lambda e: self.set(['feed', 'water_g_per_g_dry'], e.value)).classes('w-64')
                    ui.number('Moisture [% of wet feed]', value=feed.get('moisture_wt_percent'), min=0, max=99.9,
                              on_change=lambda e: self.set(['feed', 'moisture_wt_percent'], e.value)).classes('w-64')
                    ui.number('Ash [wt %, dry]', value=feed.get('ash_wt_percent'), min=0, max=100,
                              on_change=lambda e: self.set(['feed', 'ash_wt_percent'], e.value)).classes('w-64') \
                        .tooltip('Bookkeeping only; elements + ash must then sum to 100')
                    ui.select({'dry': 'dry (incl. ash)', 'dry_ash_free': 'dry, ash free (sums to 100)'},
                              label='Analysis basis', value=feed.get('basis', 'dry'),
                              on_change=lambda e: self.set(['feed', 'basis'],
                                                           None if e.value == 'dry' else e.value)).classes('w-64')
            mode.on_value_change(lambda e: self._switch_feed(e.value))
            ui.label('Added salts (they supply atoms; the products decide where the atoms go)') \
                .classes('text-body2 q-mt-sm')
            self._salts()

    def _switch_feed(self, mode):
        feed = self.data.setdefault('feed', {})
        if mode == 'species':
            feed.pop('elements_wt_percent', None)
            feed.setdefault('species_moles', {'cellulose': 1.0})
        else:
            feed.pop('species_moles', None)
            feed.setdefault('elements_wt_percent', {'C': 50.0, 'H': 6.0, 'O': 44.0})
        self.form.refresh()
        self.changed()

    def _mapping(self, path, key_label, value_label, key_options=None):
        """Rows of key + number (element wt%, species moles, atomic weights)."""
        items = list((get(self.data, path) or {}).items())

        def commit():
            self.set(path, {k: v for k, v in items if k and v is not None} or None)

        @ui.refreshable
        def rows():
            for i, (key, value) in enumerate(items):
                with ui.row().classes('items-center no-wrap gap-2'):
                    if key_options:
                        options = sorted(set(key_options) | ({key} if key else set()))
                        ui.select(options, value=key or None, with_input=True, label=key_label,
                                  on_change=lambda e, i=i: (items.__setitem__(i, (e.value, items[i][1])),
                                                            commit())).classes('w-56')
                    else:
                        ui.input(key_label, value=key,
                                 on_change=lambda e, i=i: (items.__setitem__(i, (e.value.strip(), items[i][1])),
                                                           commit())).classes('w-24')
                    ui.number(value_label, value=value,
                              on_change=lambda e, i=i: (items.__setitem__(i, (items[i][0], e.value)),
                                                        commit())).classes('w-36')
                    ui.button(icon='delete', on_click=lambda i=i: (items.pop(i), commit(), rows.refresh())) \
                        .props('flat dense round')
            ui.button(f'Add {key_label.lower()}', icon='add',
                      on_click=lambda: (items.append(('', 0.0)), rows.refresh())).props('flat dense')
        rows()

    def _salts(self):
        items = [[name, salt.get('mol_per_g_dry_feed', 0.0), format_elements(salt.get('elements'))]
                 for name, salt in (get(self.data, ['feed', 'salts']) or {}).items()]

        def commit():
            salts = {}
            for name, amount, elements in items:
                if name:
                    try:
                        salts[name] = {'mol_per_g_dry_feed': amount or 0.0, 'elements': parse_elements(elements)}
                    except ValueError:
                        return
            self.set(['feed', 'salts'], salts or None)

        def element_check(i):
            def validate(text):
                try:
                    parse_elements(text)
                except ValueError as exc:
                    return str(exc)
                items[i][2] = text
                commit()
                return None
            return validate

        @ui.refreshable
        def rows():
            for i, (name, amount, elements) in enumerate(items):
                with ui.row().classes('items-center no-wrap gap-2'):
                    ui.input('Salt', value=name,
                             on_change=lambda e, i=i: (items[i].__setitem__(0, e.value.strip()), commit())) \
                        .classes('w-28')
                    ui.number('mol / g dry feed', value=amount, format='%.3g',
                              on_change=lambda e, i=i: (items[i].__setitem__(1, e.value), commit())).classes('w-40')
                    ui.input('Elements, e.g. Na=1, Cl=1', value=elements, validation=element_check(i)) \
                        .classes('w-64')
                    ui.button(icon='delete', on_click=lambda i=i: (items.pop(i), commit(), rows.refresh())) \
                        .props('flat dense round')
            ui.button('Add salt', icon='add',
                      on_click=lambda: (items.append(['', 1e-5, '']), rows.refresh())).props('flat dense')
        rows()

    # phases
    def _phases(self):
        with self._section('Products allowed in each phase',
                           'A listed species may form; it is not forced. Ions only in the aqueous phase. '
                           'Type to search; the Species tab explains why a species is unavailable.'):
            for phase, label in PHASE_NAMES.items():
                current = list(get(self.data, ['phases', phase]) or [])
                options = list(dict.fromkeys(self.state.phase_options[phase] + current))
                ui.select(options, label=label, value=current, multiple=True, with_input=True,
                          on_change=lambda e, p=phase: self.set(['phases', p], list(e.value or []))) \
                    .classes('w-full').props('use-chips')

    # aqueous / constraints / solver
    def _aqueous(self):
        data = self.data
        with ui.card().classes('grow'):
            ui.label('Aqueous phase').classes('text-subtitle1 text-weight-medium')
            ui.select({'electrolyte': 'electrolyte (ePC-SAFT, ions)', 'molecular': 'molecular (PC-SAFT)'},
                      label='Model', value=get(data, ['aqueous', 'model'], 'electrolyte'),
                      on_change=lambda e: self.set(['aqueous', 'model'], e.value)).classes('w-full')
            ui.number('Removed at and above [°C]', value=get(data, ['aqueous', 'cutoff_C']),
                      on_change=lambda e: self.set(['aqueous', 'cutoff_C'], e.value)).classes('w-full') \
                .props('clearable').tooltip('A modelling rule, not a computed phase boundary. Empty: no cut-off')
            ui.switch('Aqueous phase present', value=get(data, ['aqueous', 'enabled'], True) is not False,
                      on_change=lambda e: self.set(['aqueous', 'enabled'], None if e.value else False))

    def _constraints(self):
        data = self.data
        with ui.card().classes('grow'):
            ui.label('Carbon constraints').classes('text-subtitle1 text-weight-medium')
            exact = get(data, ['constraints', 'gas_carbon_fraction'])
            minimum = get(data, ['constraints', 'gas_carbon_min_fraction'])
            mode = 'exact' if exact is not None else 'min' if minimum is not None else 'free'
            toggle = ui.toggle({'free': 'free', 'exact': 'exactly', 'min': 'at least'}, value=mode).props('dense')
            value = ui.number('Gas carbon fraction (0-1)', value=exact if exact is not None else minimum,
                              min=0, max=1, step=0.05).classes('w-full')

            def update(_=None):
                put(data, ['constraints', 'gas_carbon_fraction'], value.value if toggle.value == 'exact' else None)
                put(data, ['constraints', 'gas_carbon_min_fraction'], value.value if toggle.value == 'min' else None)
                self.changed()
            toggle.on_value_change(update)
            value.on_value_change(update)
            value.bind_visibility_from(toggle, 'value', backward=lambda v: v != 'free')
            ui.number('Solid carbon fraction (0-1, empty = free)',
                      value=get(data, ['constraints', 'solid_carbon_fraction']), min=0, max=1,
                      on_change=lambda e: self.set(['constraints', 'solid_carbon_fraction'], e.value)) \
                .classes('w-full').props('clearable')
            ui.select({'char': 'char only', 'all_solids': 'all solids'}, label='Solid carbon counts',
                      value=get(data, ['constraints', 'solid_carbon_target'], 'char'),
                      on_change=lambda e: self.set(['constraints', 'solid_carbon_target'],
                                                   None if e.value == 'char' else e.value)).classes('w-full')

    def _solver(self):
        data = self.data
        with ui.card().classes('grow'):
            ui.label('Solver').classes('text-subtitle1 text-weight-medium')
            ui.select(list(SOLVERS), label='Method', value=get(data, ['solver', 'method'], 'ipopt'),
                      on_change=lambda e: self.set(['solver', 'method'], e.value)).classes('w-full')
            ui.number('Multistart (extra starts)', value=get(data, ['solver', 'multistart'], 0), min=0, step=1,
                      precision=0, on_change=lambda e: self.set(['solver', 'multistart'],
                                                                int(e.value) if e.value else None)) \
                .classes('w-full').tooltip('Runs from several starting points and keeps the lowest Gibbs energy')
            ui.switch('Continuation (start from the previous point)',
                      value=get(data, ['sweep', 'continuation'], True) is not False,
                      on_change=lambda e: self.set(['sweep', 'continuation'], None if e.value else False)) \
                .tooltip('Off: every point starts fresh, and all points can run in parallel')

    # advanced
    def _advanced(self):
        data = self.data
        with ui.expansion('Advanced settings', icon='tune').classes('w-full'):
            with ui.row().classes('w-full gap-6 items-start'):
                with ui.column():
                    ui.label('Starting point').classes('text-weight-medium')
                    ui.select(list(INITIALIZATIONS), label='Initialization',
                              value=get(data, ['initial_guess', 'initialization'], 'auto'),
                              on_change=lambda e: self.set(['initial_guess', 'initialization'],
                                                           None if e.value == 'auto' else e.value)).classes('w-56')
                    self._optional_number('Initial gas water fraction', ['initial_guess', 'initial_gas_water_fraction'])
                    self._optional_number('Aqueous water floor [mol/mol C]', ['aqueous', 'water_floor'])
                with ui.column():
                    ui.label('Solver details').classes('text-weight-medium')
                    self._optional_number('Max iterations', ['solver', 'max_iterations'], integer=True)
                    self._optional_number('Random seed (multistart)', ['solver', 'seed'], integer=True)
                    ui.select({'default': 'default (with multistart)', 'on': 'on', 'off': 'off'},
                              label='Stability test', value={True: 'on', False: 'off'}.get(
                                  get(data, ['solver', 'stability_test']), 'default'),
                              on_change=lambda e: self.set(['solver', 'stability_test'],
                                                           {'on': True, 'off': False}.get(e.value))).classes('w-56')
                    self._optional_number('Energy scale [J/mol C]', ['solver', 'energy_scale'])
                    self._optional_number('slsqp-fd: ftol', ['solver', 'fd_ftol'])
                    self._optional_number('slsqp-fd: step', ['solver', 'fd_eps'])
                with ui.column():
                    ui.label('Extra temperatures (grid only)').classes('text-weight-medium')
                    refine = get(data, ['sweep', 'refine']) or {}
                    fields = {}
                    for key, label in (('start_C', 'from [°C]'), ('stop_C', 'to [°C]'), ('step_C', 'step [°C]')):
                        fields[key] = ui.number(label, value=refine.get(key)).classes('w-40').props('clearable')

                    def update_refine(_=None):
                        values = {k: f.value for k, f in fields.items()}
                        self.set(['sweep', 'refine'], values if all(v is not None for v in values.values()) else None)
                    for f in fields.values():
                        f.on_value_change(update_refine)
                with ui.column():
                    ui.label('Carbon in other phases (exact fraction)').classes('text-weight-medium')
                    for phase in ('aqueous', 'organic_liquid'):
                        self._optional_number(PHASE_NAMES[phase], ['constraints', 'phase_carbon_fractions', phase])
                    ui.label('Minimum amounts [mol/mol C]').classes('text-weight-medium q-mt-sm')
                    minimum = get(data, ['constraints', 'minimum_amounts']) or {}
                    text = '\n'.join(f'{p}, {s}, {v:g}' for p, names in minimum.items() for s, v in names.items())
                    ui.textarea('One "phase, species, value" per line', value=text,
                                validation=self._minimum_validator).classes('w-72').props('autogrow')
            with ui.row().classes('w-full gap-6 items-start'):
                with ui.column():
                    ui.label('Atomic weight overrides [g/mol]').classes('text-weight-medium')
                    ui.label('Change feed and balance arithmetic only').classes('text-caption')
                    self._mapping(['atomic_weights'], 'Element', 'g/mol')
                with ui.column():
                    ui.label('Comparison data in plots').classes('text-weight-medium')
                    files = self.state.comparison_files(self.document.base_dir())
                    current = get(data, ['comparison', 'file'])
                    options = {'': '(none)', **files, **({current: current} if current and current not in files else {})}
                    ui.select(options, label='CSV file', value=current or '',
                              on_change=lambda e: self.set(['comparison', 'file'], e.value or None)).classes('w-80')
                    ui.select(QUANTITY_CHOICES, label='Its unit',
                              value=get(data, ['comparison', 'quantity'], 'mol_per_kg_dry_feed'),
                              on_change=lambda e: self.set(['comparison', 'quantity'], e.value)).classes('w-80')
                    ui.input('Label', value=get(data, ['comparison', 'label'], ''),
                             on_change=lambda e: self.set(['comparison', 'label'], e.value)).classes('w-80')

    def _optional_number(self, label, path, integer=False):
        value = get(self.data, path)

        def change(e):
            v = e.value
            self.set(path, (int(v) if integer else v) if v is not None else None)
        return ui.number(label, value=value, on_change=change, precision=0 if integer else None) \
            .classes('w-56').props('clearable')

    def _minimum_validator(self, text):
        out = {}
        for line in (text or '').splitlines():
            if not line.strip():
                continue
            parts = [p.strip() for p in line.split(',')]
            if len(parts) != 3:
                return f'"{line.strip()}": write phase, species, value'
            phase, species, value = parts
            if phase not in PHASE_NAMES:
                return f'Unknown phase "{phase}" (use {", ".join(PHASE_NAMES)})'
            try:
                out.setdefault(phase, {})[species] = float(value)
            except ValueError:
                return f'"{value}" is not a number'
        self.set(['constraints', 'minimum_amounts'], out or None)
        return None
