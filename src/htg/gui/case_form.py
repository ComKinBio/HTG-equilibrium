"""The Case tab without NiceGUI: wording, field parsing, notices and derived numbers.

Everything the Case tab shows that is computed (point count, the feed numbers,
the carbon conflict notice, the Initial guesses tooltip, the first problem for
the status line) comes from the functions here, so it can be tested without a
browser.
"""
import math
from pathlib import Path

from .. import units
from ..case import PHASE_LABELS, Case, Feed
from ..feed_v2 import FeedError, feed_summary

# ---------------------------------------------------------------- wording (spec section 5.1)
FEED_CAPTION = 'The wet feed going into the reactor: solids plus water.'
FEED_CAPTION_SPECIES = ('The wet feed going into the reactor: species from the database plus water. '
                        'No ash in this mode.')
PRODUCTS_CAPTION = ('Species that may form in each phase. Type to search. Ions can only go in the '
                    'aqueous phase.')
CARBON_CAPTION = 'Share of the feed carbon that ends up in each phase. Free = the solver decides.'
CUTOFF_CAPTION = ('This is a modelling choice. If you know that the aqueous phase is definitely not '
                  'present above a certain temperature, then marking it here helps the equilibrium '
                  'solver. It is not necessary to input a value.')
METHOD_TOOLTIP = ('IPOPT is the default and usually the fastest. Try SLSQP if IPOPT fails at some '
                  'points. trust-constr is slower but sometimes more robust. "SLSQP, finite '
                  'differences" is the original method, kept for comparison with old results.')
CONTINUATION_TOOLTIP = ('Start each point from the solution at the previous temperature (same '
                        'pressure), or the previous point in a list. Usually faster and smoother '
                        'along a sweep. Turn off to start every point fresh; then all points can '
                        'run in parallel.')
ELEMENTS_CAPTION = 'Scaled to sum to 100 % if needed'
SPECIES_FOOTER = 'Only species with complete data are listed'
ADVANCED_CAPTION = 'Leave a field empty to use the default'
INITIAL_COMPOSITION_CAPTION = 'Automatic picks liquid below the aqueous cut-off and vapour above it.'
SEED_CAPTION = 'Set a number to get the same starting points every run.'
ENERGY_SCALE_CAPTION = 'Scales the Gibbs energy for the optimizer. Change only if the solver struggles.'
DISTINCT_CAPTION = 'Two results from different starting points count as different if either limit below is exceeded.'
COMPOSITION_TOL_CAPTION = 'Different if any species amount differs by more than this.'
GIBBS_TOL_CAPTION = 'Different if the Gibbs energies differ by more than this.'
ATOMIC_WEIGHTS_CAPTION = 'Override the standard values [g / mol]. Affects the feed and element balances only.'
MEASURED_CAPTION = ('The file needs a temperature column (T_C), then one column per species named '
                    'phase:species, e.g. gas:methane.')
CHECK_CAPTION = 'The checks need no equation of state. Passing does not guarantee that the solver converges.'

INPUT_FORMATS = [  # Format, Type, Gives
    ('List', '300, 350, 400', '300, 350, 400'),
    ('Range (start:stop:step)', '300:400:50', '300, 350, 400'),
    ('Evenly spaced (start:stop n=count)', '300:400 n=5', '300, 325, 350, 375, 400'),
    ('Mixed', '250, 300:400:50', '250, 300, 350, 400'),
]

PHASE_NAMES = dict(PHASE_LABELS)
CARBON_PHASES = ('gas', 'solid', 'aqueous', 'organic_liquid')
CARBON_CHOICES = {'free': 'free', 'exact': 'exactly', 'min': 'at least', 'max': 'at most'}
CARBON_RANGE = 'range'  # min and max together: only offered when the file has it
CARBON_RANGE_LABEL = 'between'
SOLID_CARBON = {'char': 'char (solid carbon) only', 'all_solids': 'all solids'}
AQUEOUS_MODELS = {'electrolyte': 'Electrolyte (ePC-SAFT, ions allowed)', 'molecular': 'Molecular (PC-SAFT, no ions)'}
METHODS = {'ipopt': 'IPOPT (default)', 'slsqp': 'SLSQP', 'trust-constr': 'trust-constr',
           'slsqp-fd': 'SLSQP, finite differences (original)'}
INITIAL_COMPOSITIONS = {'auto': 'Automatic', 'liquid': 'Water mostly liquid (aqueous)',
                        'vapor': 'Water mostly vapour (gas)', 'lp': 'From a linear program'}
STABILITY_TEST = {'default': 'On when using several initial guesses', 'on': 'Always on', 'off': 'Off'}
QUANTITY_CHOICES = {'mol_per_kg_dry_feed': 'mol / kg dry feed', 'mol_per_mol_C': 'mol / mol dry-feed C',
                    'g_per_kg_dry_feed': 'g / kg dry feed'}
MODE_TOGGLE = {'grid': 'Every T at every P', 'points': 'List of (T, P) points'}
FEED_TOGGLE = {'elements': 'Elemental analysis', 'species': 'Species amounts'}

PERIODIC_TABLE = (
    'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se '
    'Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb '
    'Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm '
    'Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og').split()


# ---------------------------------------------------------------- case data
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


def tidy(value, digits=12):
    """Drop floating-point noise from converted values: 299.99999999999994 -> 300."""
    return float(f'{value:.{digits}g}')


# ---------------------------------------------------------------- one number in a text field
def parse_number(text, label=None, *, required=False, minimum=None, maximum=None, above=None,
                 integer=False):
    """Text of a number field -> float (int if ``integer``), or None when empty and not required.

    minimum/maximum are inclusive limits, ``above`` an exclusive lower limit.
    Raises ValueError with a message for under the field; with ``label`` the
    message starts with the field name (for the status line).
    """
    def fail(message):
        raise ValueError(problem_text(label, message) if label else message)

    text = str(text if text is not None else '').strip()
    if not text:
        if required:
            fail('Give a value.')
        return None
    try:
        value = float(text)
    except ValueError:
        value = math.nan
    if not math.isfinite(value):
        fail(f'"{text}" is not a number.')
    if integer:
        if not value.is_integer():
            fail('Give a whole number.')
        value = int(value)
    if minimum is not None and value < minimum:
        fail(f'Give a value of {units.format_number(minimum)} or more.')
    if maximum is not None and value > maximum:
        fail(f'Give a value of {units.format_number(maximum)} or less.')
    if above is not None and value <= above:
        fail(f'Give a value above {units.format_number(above)}.')
    return value


def number_text(value, unit=None):
    """A stored number as field text ('' for None)."""
    if value is None:
        return ''
    return units.format_number(value, unit)


def problem_text(label, message):
    """'Ash', "Ash can't be ..." -> the message alone; otherwise 'Label: message'."""
    name = label.split(' [')[0].split(' · ')[-1].strip()
    if not name or message.lower().startswith(name.lower()):
        return message
    if message[:1].isupper() and message[1:2].islower():
        message = message[0].lower() + message[1:]
    return f'{name}: {message}'


# ---------------------------------------------------------------- temperature and pressure
def _values(value):
    if value is None:
        return []
    return list(value) if isinstance(value, (list, tuple)) else [value]


def condition_temperatures(conditions):
    """Temperatures [°C] typed in [conditions] (grid or points), without the cut-off."""
    conditions = conditions or {}
    if conditions.get('points') is not None:
        return [float(t) for t, _ in conditions['points']]
    return [float(t) for t in _values(conditions.get('temperature_C'))]


def point_count(conditions, cutoff_C=None):
    """Number of (T, P) points the case will solve, or None when the conditions are incomplete.

    Counted as the calculation does (a cut-off inside the temperature range adds a point).
    """
    conditions = conditions or {}
    try:
        case = Case(temperature_C=conditions.get('temperature_C'), pressure_bar=conditions.get('pressure_bar'),
                    points=conditions.get('points'), feed=Feed(), phases={}, aqueous_cutoff_C=cutoff_C)
        return len(case.grid())
    except (ValueError, TypeError, KeyError):
        return None


def plural(count, word, words=None):
    return f'{count} {word if count == 1 else (words or word + "s")}'


def point_count_text(count):
    """'13 points', '1 point'; '' when unknown."""
    return '' if count is None else plural(count, 'point')


def ready_text(points, pieces):
    return f'Ready to run: {plural(points, "point")} in {plural(pieces, "independent piece")}.'


def initial_guesses_tooltip(count):
    """The live tooltip of the Initial guesses field."""
    guesses = 'initial guess' if count == 1 else 'initial guesses'
    return (f'You will start from {count} {guesses}. The optimization algorithm finds local equilibrium '
            'points, but these may not be the global equilibrium. More starting points increase '
            'confidence in finding the global equilibrium, at a higher computational cost. If several '
            'distinct local equilibrium points are found, they are listed on the Results tab so you '
            'can compare them.')


def conflict_notice(data, temperature_unit='C'):
    """The Aqueous phase card's conflict notice, or None.

    Shown when (1) the cut-off has a value, (2) the aqueous carbon choice is
    exactly or at least with a value above 0, and (3) at least one temperature
    of the case is at or above the cut-off.
    """
    cutoff = get(data, ['aqueous', 'cutoff_C'])
    bound = get(data, ['constraints', 'carbon', 'aqueous']) or {}
    if cutoff is None or not isinstance(bound, dict):
        return None
    kind = 'exact' if 'exact' in bound else 'min' if 'min' in bound else None
    if kind is None or not isinstance(bound[kind], (int, float)) or bound[kind] <= 0:
        return None
    try:
        temperatures = condition_temperatures(data.get('conditions'))
        cutoff = float(cutoff)
    except (TypeError, ValueError):
        return None
    if not any(t >= cutoff for t in temperatures):
        return None
    choice = 'exactly' if kind == 'exact' else 'at least'
    percent = units.format_number(units.fraction_to_percent(bound[kind]))
    T = f'{units.display(cutoff, temperature_unit)} {units.label(temperature_unit)}'
    return (f'Aqueous carbon is set to {choice} {percent} % (Carbon split). At {T} and above there is no '
            f'aqueous phase, so that setting is ignored there. Below {T} it applies.')


# ---------------------------------------------------------------- carbon split
def carbon_choice(bound):
    """(choice, first fraction, second fraction) for a bound of [constraints.carbon].

    choice is a key of CARBON_CHOICES, or CARBON_RANGE for min and max together
    (first = min, second = max).
    """
    if not isinstance(bound, dict) or not bound:
        return 'free', None, None
    if 'exact' in bound:
        return 'exact', bound['exact'], None
    if 'min' in bound and 'max' in bound:
        return CARBON_RANGE, bound['min'], bound['max']
    if 'min' in bound:
        return 'min', bound['min'], None
    return 'max', bound.get('max'), None


def carbon_bound(choice, first=None, second=None):
    """The [constraints.carbon] entry for a row of the Carbon split card (None = free)."""
    if choice == 'free':
        return None
    if choice == CARBON_RANGE:
        return {'min': first, 'max': second}
    return {choice: first}


# ---------------------------------------------------------------- feed
def element_symbol(text):
    """'c' -> 'C', ' na ' -> 'Na'."""
    text = str(text or '').strip()
    return text[:1].upper() + text[1:].lower()


def element_error(symbol, existing=()):
    """Why ``symbol`` cannot be added to an element table (None when it can)."""
    if not symbol:
        return 'Type an element symbol, e.g. Cl.'
    if symbol not in PERIODIC_TABLE:
        return f'"{symbol}" is not an element symbol.'
    if symbol in existing:
        return f'{symbol} is already in the table.'
    return None


def entered_sum_notice(entered_sum):
    """The scaling notice, or None when the entered percentages add up to 100.00."""
    text = f'{entered_sum:.2f}'
    if text == '100.00':
        return None
    return (f'The element percentages add up to {text} %. They have been scaled to 100 % '
            '(the Used column).')


def feed_error_field(message):
    """The field a FeedError message is about: 'total', 'solids', 'ash', 'elements',
    'element:X', 'species', 'species:X' (or 'feed')."""
    for prefix, field in (('Total mass', 'total'), ('Solids', 'solids'), ('Ash', 'ash'),
                          ('Elements', 'elements'), ('Unknown element', 'elements')):
        if message.startswith(prefix):
            return field
    if message.startswith('Element '):
        return 'element:' + message[len('Element '):].split(':')[0]
    if message.startswith('Species amounts: water'):
        return 'species:water'
    if message.startswith('Species amounts'):
        return 'species'
    if message.startswith('Species '):
        return 'species:' + message[len('Species '):].split(':')[0]
    if message.startswith('Unknown feed species'):
        return 'species:' + message.split("'")[1] if "'" in message else 'species'
    if '(species ' in message:
        return 'species:' + message.rsplit('(species ', 1)[1].rstrip(').')
    return 'feed'


def check_feed(feed, atomic_weights=None, species_elements=None):
    """(summary, field, message): the feed summary, or None with the first problem."""
    try:
        return feed_summary(feed or {}, atomic_weights, species_elements), None, None
    except FeedError as exc:
        return None, feed_error_field(str(exc)), str(exc)


def convert_share(value, old_unit, new_unit, whole_kg):
    """A '%' or 'kg' amount in the other unit; ``whole_kg`` is what the % refers to.

    Returns the value unchanged when it cannot be converted (no value or no whole).
    """
    if value is None or old_unit == new_unit or not whole_kg:
        return value
    converted = value / 100.0 * whole_kg if new_unit == 'kg' else 100.0 * value / whole_kg
    return tidy(converted, 6)


def convert_species_amount(value, old_unit, new_unit, molar_mass):
    """A species amount in 'kg' or 'mol' in the other unit (``molar_mass`` in g/mol)."""
    if value is None or old_unit == new_unit or not molar_mass:
        return value
    converted = value * 1000.0 / molar_mass if new_unit == 'mol' else value * molar_mass / 1000.0
    return tidy(converted, 6)


def format_kg(value):
    """Feed summary masses: 8.500, 0.127; small values with 3 significant figures."""
    if value == 0:
        return '0'
    if abs(value) >= 0.01:
        return f'{value:.3f}'
    return units.format_number(value, digits=3)


def format_mol(value):
    """Feed summary amounts: 471.8, 37.96 (4 significant figures)."""
    if value is None:
        return ''
    if value == 0:
        return '0'
    if abs(value) >= 1000:
        return f'{value:.0f}'
    return units.format_number(value, digits=4) if abs(value) < 1e-3 else f'{value:#.4g}'.rstrip('.')


def format_percent(value):
    """Used column and sums: two decimals."""
    return f'{value:.2f}'


# ---------------------------------------------------------------- file bar
def file_status(path, dirty, cwd=None):
    """'cases/manure.toml · unsaved changes'."""
    if path is None:
        return 'New case · not saved yet'
    path = Path(path).resolve()
    try:
        where = str(path.relative_to(Path(cwd or Path.cwd()).resolve()))
    except ValueError:
        where = str(path)
    return f'{where} · unsaved changes' if dirty else where
