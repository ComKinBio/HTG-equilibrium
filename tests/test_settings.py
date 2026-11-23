"""App preferences file: defaults, saving, and recovery from broken files."""
import logging
import tomllib

import pytest

from htg import settings
from htg.settings import DEFAULTS, Settings


def test_defaults_when_missing(tmp_path):
    prefs = Settings.load(tmp_path)
    assert prefs.path == tmp_path / 'settings.toml'
    assert prefs.as_dict() == DEFAULTS
    assert prefs.get('temperature_unit') == 'C' and prefs.get('pressure_unit') == 'bar'
    assert prefs.get('color_theme') == 'Default' and prefs.get('line_width') == 2
    assert prefs.get('input_formats_open') is False
    assert prefs.export_defaults('one') == {'width_mm': 85, 'height_mm': 65, 'dpi': 300, 'format': 'png',
                                            'shared_legend': True}
    sizes = {layout: (prefs.export_defaults(layout)['width_mm'], prefs.export_defaults(layout)['height_mm'])
             for layout in settings.EXPORT_LAYOUTS}
    assert sizes == {'one': (85, 65), 'grid': (170, 130), 'row': (250, 65), 'present': (190, 65)}
    assert not prefs.path.exists()  # nothing written until a value changes


def test_values_survive_a_restart(tmp_path):
    prefs = Settings.load(tmp_path)
    prefs.set('pressure_unit', 'MPa')
    prefs.set('temperature_unit', 'K')
    prefs.set('line_width', 1.5)
    prefs.set('input_formats_open', True)
    prefs.set_export_defaults('row', width_mm=240, format='pdf')
    again = Settings.load(tmp_path)
    assert again.get('pressure_unit') == 'MPa' and again.get('temperature_unit') == 'K'
    assert again.get('line_width') == 1.5 and again.get('input_formats_open') is True
    assert again.export_defaults('row') == {'width_mm': 240, 'height_mm': 65, 'dpi': 300, 'format': 'pdf',
                                            'shared_legend': True}
    assert again.export_defaults('one') == DEFAULTS['export_defaults']['one']
    stored = tomllib.loads(again.path.read_text())
    assert stored['export_defaults']['row']['width_mm'] == 240


def test_environment_variable(tmp_path, monkeypatch):
    monkeypatch.setenv('HTG_CONFIG_DIR', str(tmp_path / 'env'))
    assert settings.settings_path() == tmp_path / 'env' / 'settings.toml'
    assert settings.settings_path(tmp_path) == tmp_path / 'settings.toml'
    Settings.load().set('color_theme', 'Set2')
    assert Settings.load(tmp_path / 'env').get('color_theme') == 'Set2'


def test_writes_only_on_change(tmp_path):
    prefs = Settings.load(tmp_path)
    prefs.set('pressure_unit', 'bar')
    prefs.set_export_defaults('one', width_mm=85)
    assert not prefs.path.exists()


def test_invalid_values_are_refused(tmp_path):
    prefs = Settings.load(tmp_path)
    with pytest.raises(ValueError):
        prefs.set('pressure_unit', 'psi')
    with pytest.raises(ValueError):
        prefs.set('line_width', 0)
    with pytest.raises(KeyError):
        prefs.set('solver', 'ipopt')
    with pytest.raises(KeyError):
        prefs.export_defaults('poster')
    with pytest.raises(ValueError):
        prefs.set_export_defaults('one', format='jpg')
    with pytest.raises(KeyError):
        prefs.set_export_defaults('one', colour='red')
    assert not prefs.path.exists()


def test_returned_values_are_copies(tmp_path):
    prefs = Settings.load(tmp_path)
    prefs.export_defaults('one')['width_mm'] = 1
    prefs.get('export_defaults')['one']['width_mm'] = 1
    assert prefs.export_defaults('one')['width_mm'] == 85


def test_broken_file_gives_defaults(tmp_path, caplog):
    (tmp_path / 'settings.toml').write_text('pressure_unit = "MPa\n[[[', encoding='utf-8')
    with caplog.at_level(logging.WARNING):
        prefs = Settings.load(tmp_path)
    assert prefs.as_dict() == DEFAULTS
    assert 'Could not read' in caplog.text


def test_bad_values_are_ignored_one_by_one(tmp_path, caplog):
    (tmp_path / 'settings.toml').write_text(
        'pressure_unit = "MPa"\ntemperature_unit = "F"\nline_width = "thick"\nunknown = 1\n'
        'export_defaults = { one = { width_mm = 90, dpi = -1 }, poster = { width_mm = 1 } }\n', encoding='utf-8')
    with caplog.at_level(logging.WARNING):
        prefs = Settings.load(tmp_path)
    assert prefs.get('pressure_unit') == 'MPa'
    assert prefs.get('temperature_unit') == 'C' and prefs.get('line_width') == 2
    assert prefs.export_defaults('one')['width_mm'] == 90 and prefs.export_defaults('one')['dpi'] == 300
    for text in ['temperature_unit', 'line_width', 'unknown', 'dpi', 'poster']:
        assert text in caplog.text


def test_unwritable_folder_only_warns(tmp_path, caplog):
    blocker = tmp_path / 'file'
    blocker.write_text('')
    prefs = Settings.load(blocker / 'config')
    with caplog.at_level(logging.WARNING):
        prefs.set('pressure_unit', 'atm')
    assert prefs.get('pressure_unit') == 'atm'
    assert 'Could not save' in caplog.text
