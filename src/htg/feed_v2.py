"""The feed of case-file format 2, as plain functions (no databases, no Julia).

Elemental analysis::

    [feed]
    total_kg = 10.0                          # solids + water
    solids = { value = 15.0, unit = "%" }    # "%" of the total, or "kg"
    ash = { value = 8.45, unit = "%" }       # "%" of the solids, or "kg"; omit = no ash
    elements_wt_percent = { C = 30.4, H = 4.41, O = 53.58, N = 2.66, S = 0.51 }

The element percentages are of the ash-free solids, as entered; they are scaled
to sum to 100 when used.

Species amounts::

    [feed.species]
    water = { value = 4.5, unit = "kg" }     # always present; "kg" or "mol"
    cellulose = { value = 0.5, unit = "kg" }

``feed_v2_to_model`` maps either form onto the model's ``Feed`` (no change to
the science); ``feed_summary`` gives the numbers the Case tab shows.
Species need molar masses: pass ``molar_mass(name) -> g/mol`` (see
``molar_mass_function``) or ``species_elements(name) -> {element: count}``.
"""
import copy
import math
from decimal import Decimal, InvalidOperation

from .constants import ATOMIC_WEIGHTS

ELEMENTAL_KEYS = ('total_kg', 'solids', 'ash', 'elements_wt_percent')
SPECIES_KEY = 'species'
V2_FEED_KEYS = ('total_kg', 'solids', 'ash', 'species')  # keys that only exist in format 2
SOLIDS_UNITS = ('%', 'kg')
SPECIES_UNITS = ('kg', 'mol')


class FeedError(ValueError):
    """The feed cannot be used; the message names the field."""


# ---------------------------------------------------------------- helpers
def _number(value, field):
    if isinstance(value, bool):
        raise FeedError(f'{field} must be a number, not {value!r}.')
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise FeedError(f'{field} must be a number, not {value!r}.') from None
    if not math.isfinite(number):
        raise FeedError(f'{field} must be a finite number.')
    return number


def amount(entry, field, units):
    """(value, unit) of an entry ``{value = x, unit = "u"}``."""
    if not isinstance(entry, dict) or 'value' not in entry:
        raise FeedError(f'{field}: write {{ value = ..., unit = "{units[0]}" }}.')
    unknown = set(entry) - {'value', 'unit'}
    if unknown:
        raise FeedError(f'{field}: unknown key(s) {sorted(unknown)}; use value and unit.')
    unit = entry.get('unit', units[0])
    if unit not in units:
        raise FeedError(f'{field}: the unit must be ' + ' or '.join(f'"{u}"' for u in units)
                        + f', not {unit!r}.')
    return _number(entry['value'], field), unit


def feed_mode(feed):
    """'species' or 'elements' for a format-2 [feed] table."""
    return 'species' if SPECIES_KEY in feed else 'elements'


def is_v2_feed(feed):
    return isinstance(feed, dict) and any(key in feed for key in V2_FEED_KEYS)


def water_molar_mass(atomic_weights=None):
    aw = {**ATOMIC_WEIGHTS, **(atomic_weights or {})}
    return 2.0 * aw['H'] + aw['O']


def molar_mass_function(species_elements, atomic_weights=None):
    """molar_mass(name) [g/mol] from species_elements(name) -> {element: count}."""
    aw = {**ATOMIC_WEIGHTS, **(atomic_weights or {})}

    def molar_mass(name):
        elements = species_elements(name)
        missing = [e for e in elements if e not in aw]
        if missing:
            raise FeedError(f'No atomic weight for {", ".join(missing)} (species {name}).')
        return float(sum(aw[e] * count for e, count in elements.items()))
    return molar_mass


def registry_elements(registry):
    """species_elements(name) for a species Registry (molecules, solids and ions)."""
    entries = {**registry.all_molecular, **registry.solid, **registry.ion}

    def species_elements(name):
        if name not in entries:
            raise FeedError(f'Unknown feed species {name!r}.')
        return dict(entries[name]['elements'])
    return species_elements


# ---------------------------------------------------------------- elemental analysis
def scaled_elements(elements):
    """(used, entered sum): the element percentages scaled to sum to 100."""
    if not isinstance(elements, dict) or not elements:
        raise FeedError('Elements: give at least one element.')
    entered = {}
    for element, value in elements.items():
        number = _number(value, f'Element {element}')
        if number < 0:
            raise FeedError(f'Element {element}: the percentage cannot be negative.')
        entered[str(element)] = number
    total = sum(entered.values())
    if total <= 0:
        raise FeedError('Elements: the percentages add up to zero.')
    if abs(total - 100.0) <= 1e-12 * 100.0:  # already 100: keep the numbers exactly
        return dict(entered), total
    return {e: v * 100.0 / total for e, v in entered.items()}, total


def elemental_masses(feed):
    """Masses [kg] of an elemental-analysis feed: total, solids, water, ash, ash_free."""
    if 'total_kg' not in feed:
        raise FeedError('Total mass: give total_kg (solids + water).')
    total = _number(feed['total_kg'], 'Total mass')
    if total <= 0:
        raise FeedError('Total mass must be positive.')
    if 'solids' not in feed:
        raise FeedError('Solids: give solids = { value = ..., unit = "%" or "kg" }.')
    value, unit = amount(feed['solids'], 'Solids', SOLIDS_UNITS)
    solids = value / 100.0 * total if unit == '%' else value
    if solids <= 0:
        raise FeedError('Solids must be positive.')
    if solids > total * (1 + 1e-12):
        raise FeedError("Solids can't be more than the total mass.")
    solids = min(solids, total)
    ash = 0.0
    if feed.get('ash') is not None:
        value, unit = amount(feed['ash'], 'Ash', SOLIDS_UNITS)
        if value < 0:
            raise FeedError('Ash cannot be negative.')
        ash = value / 100.0 * solids if unit == '%' else value
        if ash > solids * (1 + 1e-12):
            raise FeedError("Ash can't be more than the solids.")
    return {'total': total, 'solids': solids, 'water': total - solids, 'ash': ash,
            'ash_free': solids - ash}


def ash_percent(feed):
    """Ash as % of the solids, or None when the feed gives no ash."""
    if feed.get('ash') is None:
        return None
    value, unit = amount(feed['ash'], 'Ash', SOLIDS_UNITS)
    if unit == '%':
        return value
    return 100.0 * elemental_masses(feed)['ash'] / elemental_masses(feed)['solids']


# ---------------------------------------------------------------- species amounts
def species_masses(species, molar_mass):
    """{name: (kg, mol)} of a species-amounts feed (water included)."""
    if not isinstance(species, dict) or 'water' not in species:
        raise FeedError('Species amounts: water is always present; give water = { value = ..., unit = "kg" }.')
    out = {}
    for name, entry in species.items():
        value, unit = amount(entry, f'Species {name}', SPECIES_UNITS)
        if value < 0:
            raise FeedError(f'Species {name}: the amount cannot be negative.')
        M = molar_mass(name)
        if unit == 'kg':
            out[name] = (value, 1000.0 * value / M)
        else:
            out[name] = (value * M / 1000.0, value)
    if not any(kg > 0 for name, (kg, _) in out.items() if name != 'water'):
        raise FeedError('Species amounts: add at least one species besides water.')
    return out


def resolve_species(species, molar_mass):
    """(species_moles, water_g_per_g_dry) of a species-amounts feed."""
    masses = species_masses(species, molar_mass)
    moles = {name: mol for name, (_, mol) in masses.items() if name != 'water'}
    dry_kg = sum(kg for name, (kg, _) in masses.items() if name != 'water')
    return moles, masses['water'][0] / dry_kg


# ---------------------------------------------------------------- mapping to the model
def feed_v2_to_model(feed, salts=None, molar_mass=None):
    """The model ``Feed`` for a format-2 [feed] table (``salts``: {name: Salt or dict}).

    Elemental analysis: elements scaled to 100 on the dry, ash-free basis,
    ash_wt_percent = 100 ash / solids, water_g_per_g_dry = water / solids.
    Species amounts: species_moles and water_g_per_g_dry = water mass /
    non-water mass; without ``molar_mass`` the conversion is left to the
    calculation (``Feed.species_amounts``), which knows the molar masses.
    """
    from .case import Feed
    salts = dict(salts or feed.get('salts') or {})
    if feed_mode(feed) == 'species':
        extra = set(feed) - {SPECIES_KEY, 'salts'}
        if extra:
            raise FeedError(f'[feed.species] replaces {sorted(extra)}; remove them.')
        if molar_mass is None:
            species_masses(feed[SPECIES_KEY], lambda name: 1.0)  # checks the entries
            return Feed(species_amounts=copy.deepcopy(feed[SPECIES_KEY]), salts=salts)
        moles, water = resolve_species(feed[SPECIES_KEY], molar_mass)
        return Feed(species_moles=moles, water_g_per_g_dry=water, salts=salts)
    unknown = set(feed) - set(ELEMENTAL_KEYS) - {'salts'}
    if unknown:
        raise FeedError(f'Unknown [feed] key(s) {sorted(unknown)}; allowed: '
                        f'{list(ELEMENTAL_KEYS)} (or [feed.species]) and [feed.salts].')
    masses = elemental_masses(feed)
    if 'elements_wt_percent' not in feed:
        raise FeedError('Elements: give elements_wt_percent (mass % of the ash-free solids).')
    used, _ = scaled_elements(feed['elements_wt_percent'])
    ash = None if feed.get('ash') is None else 100.0 * masses['ash'] / masses['solids']
    if ash is not None and feed['ash'].get('unit', '%') == '%':
        ash = _number(feed['ash']['value'], 'Ash')  # exactly as written
    return Feed(elements_wt_percent=used, basis='dry_ash_free', ash_wt_percent=ash,
                moisture_wt_percent=0.0, water_g_per_g_dry=masses['water'] / masses['solids'],
                salts=salts)


# ---------------------------------------------------------------- summary for the form
def feed_summary(feed, atomic_weights=None, species_elements=None):
    """Numbers for the live Feed summary of the Case tab.

    Returns a dict with
        mode         'elements' or 'species'
        rows         [{'name', 'kg', 'mol'}]: Total, Water, Solids, Ash, Ash-free solids
                     (species mode: Total, Water, Solids (not water), each species),
                     then each element; 'mol' is None where it has no meaning
        elements     {element: {'kg', 'mol'}} (without the water's H and O)
        used         {element: % used} (elemental mode)
        entered_sum  sum of the entered element percentages (elemental mode)
    and, by name, total_kg, water_kg, water_mol, solids_kg, ash_kg, ash_free_kg.
    Raises FeedError (a ValueError) naming the first problem.
    """
    aw = {**ATOMIC_WEIGHTS, **(atomic_weights or {})}
    M_water = water_molar_mass(aw)
    if feed_mode(feed) == 'species':
        if species_elements is None:
            raise FeedError('Species amounts need the species data (species_elements).')
        molar_mass = molar_mass_function(species_elements, aw)
        masses = species_masses(feed[SPECIES_KEY], molar_mass)
        water_kg = masses['water'][0]
        dry_kg = sum(kg for name, (kg, _) in masses.items() if name != 'water')
        element_kg = {}
        for name, (kg, mol) in masses.items():
            if name == 'water':
                continue
            for element, count in species_elements(name).items():
                element_kg[element] = element_kg.get(element, 0.0) + mol * count * aw[element] / 1000.0
        elements = {e: {'kg': kg, 'mol': 1000.0 * kg / aw[e]} for e, kg in element_kg.items()}
        rows = [{'name': 'Total', 'kg': water_kg + dry_kg, 'mol': None},
                {'name': 'Water', 'kg': water_kg, 'mol': 1000.0 * water_kg / M_water},
                {'name': 'Solids (not water)', 'kg': dry_kg, 'mol': None}]
        rows += [{'name': name, 'kg': kg, 'mol': mol} for name, (kg, mol) in masses.items()
                 if name != 'water']
        rows += [{'name': e, **v} for e, v in elements.items()]
        return {'mode': 'species', 'rows': rows, 'elements': elements,
                'total_kg': water_kg + dry_kg, 'water_kg': water_kg,
                'water_mol': 1000.0 * water_kg / M_water, 'solids_kg': dry_kg,
                'ash_kg': 0.0, 'ash_free_kg': dry_kg}
    masses = elemental_masses(feed)
    used, entered_sum = scaled_elements(feed.get('elements_wt_percent') or {})
    missing = [e for e in used if e not in aw]
    if missing:
        raise FeedError(f'Unknown element(s) {", ".join(missing)}; give their atomic weights.')
    elements = {e: {'kg': masses['ash_free'] * p / 100.0,
                    'mol': 1000.0 * masses['ash_free'] * p / 100.0 / aw[e]} for e, p in used.items()}
    rows = [{'name': 'Total', 'kg': masses['total'], 'mol': None},
            {'name': 'Water', 'kg': masses['water'], 'mol': 1000.0 * masses['water'] / M_water},
            {'name': 'Solids', 'kg': masses['solids'], 'mol': None},
            {'name': 'Ash', 'kg': masses['ash'], 'mol': None},
            {'name': 'Ash-free solids', 'kg': masses['ash_free'], 'mol': None}]
    rows += [{'name': e, **v} for e, v in elements.items()]
    return {'mode': 'elements', 'rows': rows, 'elements': elements, 'used': used,
            'entered_sum': entered_sum, 'total_kg': masses['total'], 'water_kg': masses['water'],
            'water_mol': 1000.0 * masses['water'] / M_water, 'solids_kg': masses['solids'],
            'ash_kg': masses['ash'], 'ash_free_kg': masses['ash_free']}


def feed_fractions(feed_model, molar_mass=None):
    """(s, a) for a model Feed: solids / total mass and ash / solids.

    Results use them to convert per kg dry feed into per kg total feed (s) or
    per kg dry ash-free feed (a). For an elemental analysis the ash is the dry
    mass outside the elements (as in the legacy conversion). Species amounts are
    resolved with ``molar_mass(name) -> g/mol``; without it (or when a molar
    mass is unknown) s is None.
    """
    if feed_model.species_amounts is not None:
        if molar_mass is None:
            return None, 0.0
        try:
            feed_model = feed_model.resolved(molar_mass)
        except (KeyError, ValueError, TypeError, ZeroDivisionError):
            return None, 0.0
    s = 1.0 / (1.0 + feed_model.total_water_g_per_g_dry())
    if feed_model.species_moles is not None:
        return s, 0.0
    if feed_model.unmodelled_wt_percent:
        return s, (feed_model.ash_wt_percent or 0.0) / 100.0
    rest = 100.0 - sum(feed_model.dry_element_wt_percent(None).values())
    return s, max(0.0, rest) / 100.0


# ---------------------------------------------------------------- legacy feeds
UNMODELLED_MESSAGE = ("This case uses unmodelled_wt_percent, which the form can't show. "
                      "Edit the file directly.")
SPECIES_NEEDS_DATA_MESSAGE = ('This case gives the feed as species moles; converting it needs '
                              'the species data.')


def _plus_one(value):
    """1 + value, written without binary noise when value is a short decimal."""
    try:
        return float(Decimal(repr(float(value))) + 1)
    except (InvalidOperation, ValueError):
        return 1.0 + float(value)


def legacy_feed_to_v2(feed, molar_mass=None):
    """Format-2 [feed] table for a legacy [feed] table (same model numbers).

    Raises ``LegacyConversionError`` (from htg.case) when the feed cannot be
    written in format 2: unmodelled_wt_percent > 0, or a species feed without
    ``molar_mass``.

    Elemental feeds use 1 kg of solids. The ash written is 100 minus the dry
    element percentages, so that the dry elemental analysis the model uses is
    unchanged: for a 'dry' analysis without ash that sums to less than 100,
    the rest becomes ash; where the old file gave an ash value that differs
    from that rest (the old loader allowed 0.05 %), the rest is used.
    """
    from .case import Feed, LegacyConversionError
    feed = copy.deepcopy(dict(feed))
    salts = feed.pop('salts', None)
    if feed.get('unmodelled_wt_percent'):
        raise LegacyConversionError(UNMODELLED_MESSAGE)
    model = Feed(**feed)
    water = model.total_water_g_per_g_dry()
    if model.species_moles is not None:
        if molar_mass is None:
            raise LegacyConversionError(SPECIES_NEEDS_DATA_MESSAGE)
        dry_g = sum(float(n) * molar_mass(name) for name, n in model.species_moles.items())
        species = {'water': {'value': water * dry_g / 1000.0, 'unit': 'kg'}}
        species.update({name: {'value': float(n), 'unit': 'mol'}
                        for name, n in model.species_moles.items()})
        out = {'species': species}
    else:
        dry = model.dry_element_wt_percent(None)
        ash_given = model.ash_wt_percent
        divisor = 1.0
        if model.basis == 'dry' and ash_given:
            divisor = 1.0 - float(ash_given) / 100.0
        elements = dict(model.elements_wt_percent) if model.basis == 'dry_ash_free' else \
            {e: (float(v) / divisor if divisor != 1.0 else v) for e, v in model.elements_wt_percent.items()}
        rest = 100.0 - sum(dry.values())
        out = {'total_kg': _plus_one(water) if not model.moisture_wt_percent else 1.0 + water,
               'solids': {'value': 1.0, 'unit': 'kg'}}
        if ash_given is not None and abs(rest - float(ash_given)) <= 1e-9 * max(1.0, abs(rest)):
            out['ash'] = {'value': ash_given, 'unit': '%'}
        elif rest > 1e-9:
            out['ash'] = {'value': float(f'{rest:.14g}'), 'unit': '%'}
        out['elements_wt_percent'] = elements
    if salts:
        out['salts'] = salts
    return out
