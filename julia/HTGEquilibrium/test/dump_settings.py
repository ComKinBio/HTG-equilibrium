"""Write test/data/python_settings/: settings files written and read by the Python app.

    uv run python julia/HTGEquilibrium/test/dump_settings.py

Written (about 3 kB):

    written/settings.toml   what Settings.save writes after a few changes (non-default
                            units, a float line width, a theme name that needs escapes,
                            changed export defaults); Julia must read the same values and
                            write the same text for them
    written.json            Settings.load(written).as_dict() (the values above)
    mixed/settings.toml     a hand-written file: comments, another key order, invalid and
                            unknown values, an inline table
    mixed.json              {"values": what Python keeps, "warnings": the logged warnings}
"""
import json
import logging
import shutil
from pathlib import Path

from htg.settings import Settings

HERE = Path(__file__).resolve().parent
OUT = HERE / 'data' / 'python_settings'

MIXED = '''# Settings edited by hand
line_width = 1.5        # thinner lines
color_theme = "Set2"
temperature_unit = "F"
pressure_unit = "MPa"
input_formats_open = 1
solver = "ipopt"

[export_defaults]
one = { width_mm = 90, dpi = -1, format = "svg" }
poster = { width_mm = 1 }

[export_defaults.row]
height_mm = 70.5
shared_legend = "yes"
colour = "red"
'''


class Collect(logging.Handler):
    def __init__(self):
        super().__init__(logging.WARNING)
        self.messages = []

    def emit(self, record):
        self.messages.append(record.getMessage())


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    (OUT / 'mixed').mkdir(parents=True)

    prefs = Settings.load(OUT / 'written')
    prefs.set('pressure_unit', 'MPa')
    prefs.set('temperature_unit', 'K')
    prefs.set('line_width', 1.5)
    prefs.set('input_formats_open', True)
    prefs.set('color_theme', 'Set2 "bright" \\ é\ttab')
    prefs.set_export_defaults('row', width_mm=240, format='pdf', shared_legend=False)
    prefs.set_export_defaults('one', height_mm=62.5, dpi=600)
    (OUT / 'written.json').write_text(json.dumps(Settings.load(OUT / 'written').as_dict(), indent=1) + '\n',
                                      encoding='utf-8')

    (OUT / 'mixed' / 'settings.toml').write_text(MIXED, encoding='utf-8')
    handler = Collect()
    log = logging.getLogger('htg.settings')
    log.addHandler(handler)
    try:
        mixed = Settings.load(OUT / 'mixed')
    finally:
        log.removeHandler(handler)
    folder = str(OUT / 'mixed' / 'settings.toml')
    warnings = [m.replace(folder, '<path>') for m in handler.messages]
    (OUT / 'mixed.json').write_text(json.dumps({'values': mixed.as_dict(), 'warnings': warnings}, indent=1) + '\n',
                                    encoding='utf-8')
    print('Wrote', OUT)


if __name__ == '__main__':
    main()
