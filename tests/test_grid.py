"""Temperature-pressure grids, explicit points and continuation chains (no solving)."""
from pathlib import Path

import pytest

import htg

TEMPLATE = Path(htg.__file__).parent / 'templates' / 'case_template.toml'


def case(**changes):
    return htg.load_case(TEMPLATE).replace(**changes)


def test_full_grid_order():
    c = case(temperature_C=[400, 300], pressure_bar=[250, 240], aqueous_cutoff_C=None)
    assert c.grid() == [(300.0, 240.0), (400.0, 240.0), (300.0, 250.0), (400.0, 250.0)]
    assert c.temperatures() == [300.0, 400.0] and c.pressures() == [240.0, 250.0]


def test_chains_follow_temperature_per_pressure():
    c = case(temperature_C=[300, 400], pressure_bar=[240, 250], aqueous_cutoff_C=None)
    assert c.chains() == [[(300.0, 240.0), (400.0, 240.0)], [(300.0, 250.0), (400.0, 250.0)]]
    assert c.replace(continuation=False).chains() == [[p] for p in c.grid()]


def test_explicit_points_keep_order():
    c = case(temperature_C=None, pressure_bar=None, points=[[450, 300], [300, 240], [450, 300]])
    assert c.grid() == [(450.0, 300.0), (300.0, 240.0)]
    assert c.chains() == [[(450.0, 300.0), (300.0, 240.0)]]


def test_cutoff_and_refine_added_to_temperatures():
    c = case(temperature_C=[300, 400], pressure_bar=240, aqueous_cutoff_C=385,
             refine={'start_C': 350, 'stop_C': 370, 'step_C': 10})
    assert c.temperatures() == [300.0, 350.0, 360.0, 370.0, 385.0, 400.0]


def test_at_needs_pressure_when_several():
    c = case(temperature_C=[300], pressure_bar=[240, 250])
    with pytest.raises(ValueError, match='several pressures'):
        c.at(300)
    point = c.at(300, 250)
    assert (point.temperature_C, point.pressure_bar, point.points) == (300.0, 250.0, None)


def test_invalid_points():
    with pytest.raises(ValueError, match='temperature_C, pressure_bar'):
        case(temperature_C=None, pressure_bar=None, points=[[300]]).validate()
    with pytest.raises(ValueError, match='refine'):
        case(temperature_C=None, pressure_bar=None, points=[[300, 240]],
             refine={'start_C': 1, 'stop_C': 2, 'step_C': 1}).validate()


def test_toml_round_trip_with_grid_and_points(tmp_path):
    for c in (case(pressure_bar=[240.0, 300.0]),
              case(temperature_C=None, pressure_bar=None, points=[[300.0, 240.0], [400.0, 250.0]])):
        c.save(tmp_path / 'c.toml')
        assert htg.load_case(tmp_path / 'c.toml') == c


def test_points_exclude_temperature_in_file(tmp_path):
    text = TEMPLATE.read_text().replace('pressure_bar = 240.0', 'pressure_bar = 240.0\npoints = [[300, 240]]')
    (tmp_path / 'bad.toml').write_text(text)
    with pytest.raises(ValueError, match='leave out'):
        htg.load_case(tmp_path / 'bad.toml')
