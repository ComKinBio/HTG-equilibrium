"""Display units: conversions to and from the file units (°C, bar), and number formatting."""
import numpy as np
import pytest

from htg import units

# Temperatures stay well above 0 K: K <-> °C is an offset, so the relative error near 0 K grows.
VALUES = {'temperature': [-200.0, -40.0, 0.0, 0.1, 25.0, 300.0, 374.15, 573.15, 1234.5678],
          'pressure': [0.0, 1e-3, 0.1, 1.0, 1.01325, 240.0, 1234.5678, 2.4e7]}


def close(a, b, rel=1e-12):
    return abs(a - b) <= rel * max(abs(a), abs(b)) or a == b


@pytest.mark.parametrize('unit', units.UNITS)
def test_round_trips(unit):
    for value in VALUES[units.quantity(unit)]:
        assert close(units.from_file(units.to_file(value, unit), unit), value)
        assert close(units.to_file(units.from_file(value, unit), unit), value)


def test_conversions():
    assert units.to_file(573.15, 'K') == pytest.approx(300.0, rel=1e-14)
    assert units.from_file(300, 'K') == pytest.approx(573.15, rel=1e-14)
    assert units.to_file(24, 'MPa') == 240
    assert units.to_file(2.4e7, 'Pa') == 240
    assert units.from_file(240, 'Pa') == 2.4e7
    assert units.to_file(1, 'atm') == 1.01325
    assert units.from_file(1.01325, 'atm') == 1
    assert units.to_file(300, 'C') == 300 and units.to_file(300, '°C') == 300
    assert units.to_file(240, 'bar') == 240
    # between display units, through bar
    assert close(units.pressure_from_bar(units.pressure_to_bar(1, 'atm'), 'Pa'), 101325)
    assert close(units.pressure_from_bar(units.pressure_to_bar(1, 'MPa'), 'atm'), 10 / 1.01325)
    np.testing.assert_allclose(units.to_file(np.array([573.15, 673.15]), 'K'), [300, 400], rtol=1e-14)


def test_unknown_or_wrong_units():
    with pytest.raises(ValueError, match='Unknown unit'):
        units.to_file(1, 'psi')
    with pytest.raises(ValueError):
        units.temperature_to_C(1, 'bar')
    with pytest.raises(ValueError):
        units.pressure_to_bar(1, 'K')


def test_quantities_and_labels():
    assert [units.quantity(u) for u in units.UNITS] == ['temperature'] * 2 + ['pressure'] * 4
    assert units.label('C') == '°C' and units.label('MPa') == 'MPa'
    assert units.field_label('Temperatures', 'K') == 'Temperatures [K]'
    assert units.field_label('Pressures', 'bar') == 'Pressures [bar]'
    assert units.field_label('Temperatures', 'C') == 'Temperatures [°C]'
    assert units.FILE_UNITS == {'temperature': 'C', 'pressure': 'bar'}


def test_percent():
    assert units.percent_to_fraction(90) == 0.9
    assert units.fraction_to_percent(0.05) == 5
    for percent in [0, 0.5, 3, 33.3333, 90, 100]:
        assert close(units.fraction_to_percent(units.percent_to_fraction(percent)), percent)


def test_format_number():
    f = units.format_number
    assert f(300.0) == '300'
    assert f(0.1 + 0.2) == '0.3'
    assert f(573.15) == '573.15'
    assert f(299.99999999999994) == '300'
    assert f(1234.5678) == '1234.57'
    assert f(1234567) == '1.23457e6'
    assert f(1.5e-5) == '1.5e-5'
    assert f(-40) == '-40'
    assert f(0) == '0' and f(-0.0) == '0'
    assert f(2.4e7, 'Pa') == '2.4e7'
    assert f(101325, 'Pa') == '1.01325e5'
    assert f(1e5, 'Pa') == '1e5'
    assert f(24, 'MPa') == '24'
    assert units.format_list([300, 325.5, 350]) == '300, 325.5, 350'


def test_display_from_file_units():
    assert units.display(240, 'Pa') == '2.4e7'
    assert units.display(240, 'MPa') == '24'
    assert units.display(240, 'atm') == '236.862'
    assert units.display(300, 'K') == '573.15'
    assert units.display(300, 'C') == '300'
