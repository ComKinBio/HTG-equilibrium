"""App preferences: display choices only, never anything that changes a calculation.

Stored as TOML in ``settings.toml`` in the config folder, which is, in order:
the ``config_dir`` argument of ``Settings.load``, the ``HTG_CONFIG_DIR``
environment variable, or
``platformdirs.user_config_dir("htg-equilibrium")``. Example file::

    temperature_unit = "C"        # units.TEMPERATURE_UNITS
    pressure_unit = "bar"         # units.PRESSURE_UNITS
    color_theme = "Default"
    line_width = 2
    input_formats_open = false

    [export_defaults.one]         # Save plot dialog, one entry per layout
    width_mm = 85
    height_mm = 65
    dpi = 300
    format = "png"                # png, svg or pdf
    shared_legend = true

Layouts: ``one`` (one panel), ``grid`` (all four, 2 x 2), ``row`` (all four in
a row), ``present`` (phases present, in a row).

Use::

    settings = Settings.load()            # once, at start-up
    settings.get('pressure_unit')         # 'bar'
    settings.set('pressure_unit', 'MPa')  # validates and writes the file
    settings.export_defaults('row')       # {'width_mm': 250, ...} (a copy)
    settings.set_export_defaults('row', width_mm=240)

A missing, unreadable or invalid file (or a single invalid value) falls back to
the defaults with a logged warning. Failing to write is also only a warning.
"""
import copy
import logging
import os
import tomllib
from pathlib import Path

import platformdirs
import tomlkit

from . import units

log = logging.getLogger(__name__)

APP_NAME = 'htg-equilibrium'
FILE_NAME = 'settings.toml'
ENV_VAR = 'HTG_CONFIG_DIR'

EXPORT_LAYOUTS = ('one', 'grid', 'row', 'present')
EXPORT_FORMATS = ('png', 'svg', 'pdf')


def _export(width_mm, height_mm):
    return {'width_mm': width_mm, 'height_mm': height_mm, 'dpi': 300, 'format': 'png', 'shared_legend': True}


DEFAULTS = {
    'temperature_unit': 'C',
    'pressure_unit': 'bar',
    'color_theme': 'Default',
    'line_width': 2,
    'input_formats_open': False,
    'export_defaults': {'one': _export(85, 65), 'grid': _export(170, 130),
                        'row': _export(250, 65), 'present': _export(190, 65)},
}


def _positive_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and value > 0


def _positive_int(value):
    return isinstance(value, int) and not isinstance(value, bool) and value > 0


_CHECKS = {
    'temperature_unit': lambda v: v in units.TEMPERATURE_UNITS,
    'pressure_unit': lambda v: v in units.PRESSURE_UNITS,
    'color_theme': lambda v: isinstance(v, str) and v != '',
    'line_width': _positive_number,
    'input_formats_open': lambda v: isinstance(v, bool),
}
_EXPORT_CHECKS = {
    'width_mm': _positive_number,
    'height_mm': _positive_number,
    'dpi': _positive_int,
    'format': lambda v: v in EXPORT_FORMATS,
    'shared_legend': lambda v: isinstance(v, bool),
}


def settings_path(config_dir=None):
    """Path of settings.toml: in ``config_dir``, $HTG_CONFIG_DIR or the user config folder."""
    if config_dir is not None:
        return Path(config_dir) / FILE_NAME
    if os.environ.get(ENV_VAR):
        return Path(os.environ[ENV_VAR]) / FILE_NAME
    return Path(platformdirs.user_config_dir(APP_NAME, appauthor=False)) / FILE_NAME


class Settings:
    """Preferences read from and written to one settings file."""

    def __init__(self, path, values=None):
        self.path = Path(path)
        self._values = copy.deepcopy(DEFAULTS) if values is None else values

    @classmethod
    def load(cls, config_dir=None):
        """Read ``settings.toml`` from the config folder; defaults for anything missing or invalid."""
        path = settings_path(config_dir)
        values = copy.deepcopy(DEFAULTS)
        if not path.exists():
            return cls(path, values)
        try:
            with open(path, 'rb') as f:
                stored = tomllib.load(f)
        except (OSError, tomllib.TOMLDecodeError, UnicodeDecodeError) as exc:
            log.warning('Could not read %s (%s); using the default settings.', path, exc)
            return cls(path, values)
        for key, value in stored.items():
            if key == 'export_defaults':
                _merge_export(values['export_defaults'], value, path)
            elif key not in _CHECKS:
                log.warning('Ignoring unknown setting "%s" in %s.', key, path)
            elif _CHECKS[key](value):
                values[key] = value
            else:
                log.warning('Ignoring invalid value %r for "%s" in %s.', value, key, path)
        return cls(path, values)

    def get(self, key):
        """A preference value (export_defaults: use export_defaults(layout))."""
        if key not in DEFAULTS:
            raise KeyError(f'Unknown setting "{key}"')
        return copy.deepcopy(self._values[key])

    def set(self, key, value):
        """Change a preference and write the file if the value changed."""
        if key not in _CHECKS:
            raise KeyError(f'Unknown setting "{key}"' +
                           (' (use set_export_defaults)' if key == 'export_defaults' else ''))
        if not _CHECKS[key](value):
            raise ValueError(f'Invalid value {value!r} for setting "{key}"')
        if self._values[key] != value:
            self._values[key] = value
            self.save()

    def export_defaults(self, layout):
        """{'width_mm', 'height_mm', 'dpi', 'format', 'shared_legend'} for a layout (a copy)."""
        if layout not in EXPORT_LAYOUTS:
            raise KeyError(f'Unknown export layout "{layout}"; use one of {", ".join(EXPORT_LAYOUTS)}')
        return dict(self._values['export_defaults'][layout])

    def set_export_defaults(self, layout, **values):
        """Remember Save plot values for a layout, e.g. set_export_defaults('one', width_mm=90)."""
        current = self.export_defaults(layout)
        for key, value in values.items():
            if key not in _EXPORT_CHECKS:
                raise KeyError(f'Unknown export setting "{key}"')
            if not _EXPORT_CHECKS[key](value):
                raise ValueError(f'Invalid value {value!r} for export setting "{key}"')
        if any(current[k] != v for k, v in values.items()):
            self._values['export_defaults'][layout].update(values)
            self.save()

    def as_dict(self):
        return copy.deepcopy(self._values)

    def save(self):
        """Write the file (all values). Errors are logged, not raised."""
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            temporary = self.path.with_suffix('.toml.tmp')
            temporary.write_text(tomlkit.dumps(self._values), encoding='utf-8')
            temporary.replace(self.path)
        except OSError as exc:
            log.warning('Could not save the settings to %s (%s).', self.path, exc)


def _merge_export(target, stored, path):
    if not isinstance(stored, dict):
        log.warning('Ignoring invalid "export_defaults" in %s.', path)
        return
    for layout, values in stored.items():
        if layout not in target or not isinstance(values, dict):
            log.warning('Ignoring export_defaults.%s in %s.', layout, path)
            continue
        for key, value in values.items():
            if key in _EXPORT_CHECKS and _EXPORT_CHECKS[key](value):
                target[layout][key] = value
            else:
                log.warning('Ignoring export_defaults.%s.%s = %r in %s.', layout, key, value, path)
