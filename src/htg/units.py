"""Display units for temperature, pressure and carbon bounds.

Case files always store temperature in °C, pressure in bar (absolute) and carbon
bounds as fractions 0-1 (the *file units*). The interface shows and accepts other
units; every conversion goes through this module.

Units are named by short keys, the same strings the settings file uses:

    temperature: 'C' (default), 'K'
    pressure:    'bar' (default), 'MPa', 'Pa', 'atm'

``'°C'`` is accepted as another spelling of ``'C'``.

    to_file(573.15, 'K')        -> 300.0        (display unit -> °C or bar)
    from_file(240, 'MPa')       -> 24.0         (°C or bar -> display unit)
    format_number(2.4e7, 'Pa')  -> '2.4e7'
    display(240, 'Pa')          -> '2.4e7'      (from_file, then format_number)
    field_label('Temperatures', 'K') -> 'Temperatures [K]'
    percent_to_fraction(90)     -> 0.9

The conversions work on floats and on NumPy arrays.
"""
import math

TEMPERATURE_UNITS = ('C', 'K')
PRESSURE_UNITS = ('bar', 'MPa', 'Pa', 'atm')
UNITS = TEMPERATURE_UNITS + PRESSURE_UNITS
FILE_UNITS = {'temperature': 'C', 'pressure': 'bar'}
LABELS = {'C': '°C', 'K': 'K', 'bar': 'bar', 'MPa': 'MPa', 'Pa': 'Pa', 'atm': 'atm'}

ZERO_CELSIUS_K = 273.15
BAR_PER_ATM = 1.01325
BAR_PER_MPA = 10.0
PA_PER_BAR = 1e5
_SCIENTIFIC = {'Pa'}
SIGNIFICANT_FIGURES = 6


def unit_key(unit):
    """The canonical key of ``unit`` ('°C' -> 'C'); ValueError for an unknown unit."""
    key = 'C' if unit == '°C' else unit
    if key not in UNITS:
        raise ValueError(f'Unknown unit "{unit}". Use one of: {", ".join(LABELS.values())}.')
    return key


def quantity(unit):
    """'temperature' or 'pressure'."""
    return 'temperature' if unit_key(unit) in TEMPERATURE_UNITS else 'pressure'


def label(unit):
    """The unit as shown: 'C' -> '°C', 'MPa' -> 'MPa'."""
    return LABELS[unit_key(unit)]


def field_label(name, unit):
    """'Temperatures', 'K' -> 'Temperatures [K]'."""
    return f'{name} [{label(unit)}]'


# ---------------------------------------------------------------- conversions
def temperature_to_C(value, unit):
    unit = unit_key(unit)
    if unit == 'C':
        return value
    if unit == 'K':
        return value - ZERO_CELSIUS_K
    raise ValueError(f'"{unit}" is not a temperature unit.')


def temperature_from_C(value, unit):
    unit = unit_key(unit)
    if unit == 'C':
        return value
    if unit == 'K':
        return value + ZERO_CELSIUS_K
    raise ValueError(f'"{unit}" is not a temperature unit.')


def pressure_to_bar(value, unit):
    unit = unit_key(unit)
    if unit == 'bar':
        return value
    if unit == 'MPa':
        return value * BAR_PER_MPA
    if unit == 'Pa':
        return value / PA_PER_BAR
    if unit == 'atm':
        return value * BAR_PER_ATM
    raise ValueError(f'"{unit}" is not a pressure unit.')


def pressure_from_bar(value, unit):
    unit = unit_key(unit)
    if unit == 'bar':
        return value
    if unit == 'MPa':
        return value / BAR_PER_MPA
    if unit == 'Pa':
        return value * PA_PER_BAR
    if unit == 'atm':
        return value / BAR_PER_ATM
    raise ValueError(f'"{unit}" is not a pressure unit.')


def to_file(value, unit):
    """A value in ``unit`` -> °C (temperature units) or bar (pressure units)."""
    if quantity(unit) == 'temperature':
        return temperature_to_C(value, unit)
    return pressure_to_bar(value, unit)


def from_file(value, unit):
    """A value in °C or bar -> ``unit`` (°C for temperature units, bar for pressure units)."""
    if quantity(unit) == 'temperature':
        return temperature_from_C(value, unit)
    return pressure_from_bar(value, unit)


def percent_to_fraction(percent):
    """Carbon bounds: 90 (%) -> 0.9."""
    return percent / 100.0


def fraction_to_percent(fraction):
    """Carbon bounds: 0.9 -> 90 (%)."""
    return fraction * 100.0


# ---------------------------------------------------------------- display
def _strip(mantissa):
    return mantissa.rstrip('0').rstrip('.') if '.' in mantissa else mantissa


def format_number(value, unit=None, digits=SIGNIFICANT_FIGURES):
    """A number already in ``unit`` as text: up to ``digits`` significant figures,
    no trailing zeros. Pa (and very large or small numbers) in scientific notation
    without '+' or leading zeros: 2.4e7, 1.5e-5.

    300.0 -> '300', 0.1 + 0.2 -> '0.3', 1234567 -> '1.23457e6'.
    """
    value = float(value)
    if not math.isfinite(value):
        return str(value)
    if value == 0:
        return '0'
    scientific = unit is not None and unit_key(unit) in _SCIENTIFIC
    text = f'{value:.{digits - 1}e}' if scientific else f'{value:.{digits}g}'
    if 'e' not in text:
        return _strip(text)
    mantissa, exponent = text.split('e')
    return f'{_strip(mantissa)}e{int(exponent)}'


def format_list(values, unit=None, digits=SIGNIFICANT_FIGURES):
    """[300, 350] -> '300, 350' (numbers already in ``unit``)."""
    return ', '.join(format_number(v, unit, digits) for v in values)


def display(value, unit, digits=SIGNIFICANT_FIGURES):
    """A value in file units (°C or bar) as text in ``unit``: display(240, 'MPa') -> '24'."""
    return format_number(from_file(value, unit), unit, digits)
