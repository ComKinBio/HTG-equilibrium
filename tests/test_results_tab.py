"""Results tab (WP5) with NiceGUI's simulated user, on hand-built result tables (no Julia)."""
import asyncio
import io
import struct
import zipfile
from types import SimpleNamespace

import pandas as pd
import pytest
from nicegui import ui
from test_plot_core import make_tables

from htg import plotspec as ps
from htg import quantities as q
from htg.gui.results_tab import ResultsTab, run_feed_fractions
from htg.gui.runs import RunTables
from htg.settings import Settings

V2_CASE = '''name = "grid"
format = 2
[feed]
total_kg = 10.0
solids = { value = 15.0, unit = "%" }
ash = { value = 8.0, unit = "%" }
elements_wt_percent = { C = 50, H = 6, O = 44 }
'''


def grid_tables(folder, case_text=V2_CASE, label='grid'):
    """Two pressures (240, 300 bar) x 300, 400, 450 (failed), 500 °C; no organic liquid."""
    parts = [make_tables(P=240.0), make_tables(P=300.0)]
    joined = {k: pd.concat([getattr(p, k) for p in parts], ignore_index=True) for k in ('summary', 'species', 'phases')}
    folder.mkdir(parents=True, exist_ok=True)
    if case_text is not None:
        (folder / 'case.toml').write_text(case_text, encoding='utf-8')
    return RunTables(label, None, joined['summary'], joined['species'], joined['phases'], None, folder)


class StubState:
    """What the Results tab needs from htg.gui.app.State."""

    def __init__(self, runs, settings):
        self.runs = runs
        self.settings = settings

    def sources(self):
        return {key: tables.label for key, tables in self.runs.items()}

    def run(self, source):
        return None

    def tables(self, source):
        return self.runs.get(source)


def simulate(state, scenario):
    from nicegui.testing.user_simulation import user_simulation

    holder = {}

    def page():
        holder['tab'] = ResultsTab(state)

    async def main():
        async with user_simulation(page) as user:
            await user.open('/')
            await scenario(user, holder['tab'])

    asyncio.run(main())


@pytest.fixture
def state(tmp_path, config_dir):
    runs = {'saved:grid': grid_tables(tmp_path / 'grid'),
            'saved:old': grid_tables(tmp_path / 'old', case_text=None, label='old')}
    return StubState(runs, Settings.load(config_dir))


def test_feed_fractions_of_a_run(tmp_path):
    assert run_feed_fractions(grid_tables(tmp_path / 'a')) == pytest.approx((0.15, 0.08))
    assert run_feed_fractions(grid_tables(tmp_path / 'b', case_text=None)) is None
    legacy = 'name = "x"\n[feed]\nelements_wt_percent = { C = 50 }\nwater_g_per_g_dry = 4.0\n'
    try:
        from htg.case import legacy_to_v2  # noqa: F401  (WP2)
    except ImportError:
        assert run_feed_fractions(grid_tables(tmp_path / 'c', case_text=legacy)) is None
    else:
        assert run_feed_fractions(grid_tables(tmp_path / 'c', case_text=legacy)) == pytest.approx((0.2, 0.0))


def _figure(tab):
    figure = tab.plot.figure
    return figure.to_dict() if hasattr(figure, 'to_dict') else figure


def test_units_convert_every_value(state):
    state.settings.set('pressure_unit', 'MPa')  # the Results tab starts from the Case-tab units

    async def scenario(user, tab):
        await user.should_see('Choose a run.')
        with user:
            tab.select('saved:grid')
        await user.should_see('6 of 8 points solved')
        assert tab.spec.title == 'Amount in each phase · mol / kg dry feed · 24 MPa'
        assert list(tab.fixed.options.values()) == ['24 MPa', '30 MPa']
        layout = _figure(tab)['layout']
        assert layout['xaxis']['title']['text'] == 'Temperature [°C]'
        gas = tab.spec.panels[0].series[0]
        assert gas.x == [300.0, 400.0, 450.0, 500.0]

        with user:
            tab.x_unit.set_value('K')
        layout = _figure(tab)['layout']
        assert layout['xaxis']['title']['text'] == 'Temperature [K]'
        assert tab.spec.panels[0].series[0].x == pytest.approx([573.15, 673.15, 723.15, 773.15])
        assert _figure(tab)['data'][0]['x'][0] == pytest.approx(573.15)  # the ticks follow the data

        with user:
            tab.fixed_unit.set_value('Pa')
        assert list(tab.fixed.options.values()) == ['2.4e7 Pa', '3e7 Pa']
        assert tab.spec.title.endswith(' · 2.4e7 Pa')

        with user:
            tab.x_axis.set_value('P')
        assert tab.fixed.label == 'At temperature'
        assert list(tab.fixed.options.values()) == ['573.15 K', '673.15 K', '723.15 K', '773.15 K']
        assert tab.x_unit.value == 'Pa' and tab.fixed_unit.value == 'K'
        assert _figure(tab)['layout']['xaxis']['title']['text'] == 'Pressure [Pa]'
        assert tab.spec.panels[0].series[0].x == pytest.approx([2.4e7, 3e7])
        assert tab.spec.title.endswith(' · 573.15 K')

        # the T–P map: its tick labels are the run's values in the current units
        with user:
            tab.kind.set_value('map')
        heat = tab.spec.panels[0].heatmap
        assert heat.x == ['573.15', '673.15', '723.15', '773.15'] and heat.y == ['2.4e7', '3e7']
        assert _figure(tab)['data'][0]['x'] == heat.x
        assert state.settings.get('temperature_unit') == 'C'  # Results units are not preferences

    simulate(state, scenario)


def test_bases_and_species(state):
    async def scenario(user, tab):
        with user:
            tab.select('saved:grid')
        assert 'disable' not in tab.basis.props['options'][0]
        with user:
            tab.basis.set_value('total')
        assert tab.spec.title.startswith('Amount in each phase · mol / kg total feed')
        gas = {s.name: s for s in tab.spec.panels[0].series}
        assert gas['carbon dioxide'].y[0] == pytest.approx(10.0 * 0.15)
        assert [p.empty for p in tab.spec.panels] == [False, False, False, True]

        # species list: largest 8 ticked, water off; unticking removes the line
        assert tab.species_checks['carbon dioxide'].value and not tab.water.value
        user.find(kind=ui.checkbox, content='methane').click()
        assert 'methane' not in [s.name for s in tab.spec.panels[0].series]
        user.find(kind=ui.checkbox, content='water').click()
        assert 'water' in [s.name for s in tab.spec.panels[0].series]
        user.find('Largest 8').click()
        assert 'methane' in [s.name for s in tab.spec.panels[0].series]
        user.find('Filter').type('hydro')
        assert list(tab.species_checks) == ['hydrogen']

        # a run without a total mass: the basis falls back and is disabled with a tooltip
        with user:
            tab.select('saved:old')
        assert tab.basis.value == 'dry'
        assert tab.basis.props['options'][0]['disable'] is True
        assert tab.basis.reason.text == q.NO_TOTAL_MASS and tab.basis.reason.visible

        # gas composition options are enabled only for that plot
        assert not tab.composition.enabled
        user.find('Gas composition').click()
        assert tab.kind.value == 'gas' and tab.composition.enabled and not tab.amount.enabled

    simulate(state, scenario)


def test_colours_and_line_are_preferences(state):
    async def scenario(user, tab):
        with user:
            tab.select('saved:grid')
            tab.theme.set_value('Set2')
            tab.line.set_value(3)
        assert state.settings.get('color_theme') == 'Set2' and state.settings.get('line_width') == 3
        assert len(tab.strip.default_slot.children) == len(ps.PALETTES['Set2'])
        trace = next(t for t in _figure(tab)['data'] if t.get('name') == 'carbon dioxide')
        assert trace['line']['color'] in ps.PALETTES['Set2'] and trace['line']['width'] == 3
        assert Settings.load(state.settings.path.parent).get('color_theme') == 'Set2'

    simulate(state, scenario)


def _png_size(data):
    return struct.unpack('>II', data[16:24])


def test_save_plot_dialog(state):
    async def scenario(user, tab):
        with user:
            tab.select('saved:grid')
        user.find('Save plot…').click()
        dialog = tab.save_dialog
        assert list(dialog.layout.options.values()) == ['One panel', 'All four, 2 × 2', 'All four, in a row',
                                                        'Phases present, in a row']
        assert dialog.layout.value == 'grid' and dialog.width.value == 170 and dialog.height.value == 130
        assert dialog.caption.text == 'As shown on the Results tab'
        with user:
            dialog.layout.set_value('present')
        assert dialog.caption.text == 'Gas, Aqueous, Solid (no organic liquid at any point)'
        assert (dialog.width.value, dialog.height.value) == (190, 65)
        assert dialog.summary.text == '3 panels · 190 × 65 mm · PNG'
        assert dialog.preview.source.startswith('data:image/png;base64,')
        with user:
            dialog.width.set_value(100)
            dialog.height.set_value(50)
            dialog.dpi.set_value(100)
        assert dialog.summary.text == '3 panels · 100 × 50 mm · PNG'
        user.find('Save…').click()
        response = await user.download.next(timeout=10)
        w, h = _png_size(response.content)
        assert abs(w - 100 / 25.4 * 100) <= 1 and abs(h - 50 / 25.4 * 100) <= 1
        remembered = state.settings.export_defaults('present')
        assert remembered['width_mm'] == 100 and remembered['dpi'] == 100
        assert state.settings.export_defaults('grid')['width_mm'] == 170

        # a vector format and one panel
        user.find('Save plot…').click()
        dialog = tab.save_dialog
        with user:
            dialog.layout.set_value('one')
            dialog.format.set_value('PDF')
        assert not dialog.dpi.visible and dialog.panel.visible
        assert dialog.summary.text == '1 panel · 85 × 65 mm · PDF'
        user.find('Save…').click()
        response = await user.download.next(timeout=10)
        assert response.content.startswith(b'%PDF')

        # single-panel plots only offer one panel
        user.find('Carbon gasification efficiency').click()
        user.find('Save plot…').click()
        assert list(tab.save_dialog.layout.options) == ['one']
        assert not tab.save_dialog.panel.visible

    simulate(state, scenario)


def test_tables_and_download(state):
    async def scenario(user, tab):
        with user:
            tab.select('saved:grid')
        assert len(tab.summary_table.rows) == 8
        user.find('Download all as CSV').click()
        response = await user.download.next(timeout=5)
        with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
            assert sorted(archive.namelist()) == ['phases.csv', 'species.csv', 'summary.csv']
            summary = pd.read_csv(archive.open('summary.csv'))
        assert len(summary) == 8

    simulate(state, scenario)


def test_map_needs_a_grid(tmp_path, config_dir):
    one = make_tables()
    tables = RunTables('line', None, one.summary, one.species, one.phases, None, tmp_path)
    state = StubState({'saved:line': tables}, Settings.load(config_dir))

    async def scenario(user, tab):
        with user:
            tab.select('saved:line')
        assert tab.kind.props['options'][-1]['disable'] is True
        assert 'needs at least 2 temperatures' in tab.kind.reason.text

    simulate(state, scenario)


def test_no_points_yet(tmp_path, config_dir):
    empty = SimpleNamespace(summary=pd.DataFrame(columns=['T_C', 'P_bar', 'status']), species=pd.DataFrame(),
                            phases=pd.DataFrame())
    tables = RunTables('new', None, empty.summary, empty.species, empty.phases, None, tmp_path)
    state = StubState({'live:1': tables}, Settings.load(config_dir))

    async def scenario(user, tab):
        with user:
            tab.select('live:1')
        await user.should_see('No points solved yet.')
        assert not tab.save_button.enabled

    simulate(state, scenario)
