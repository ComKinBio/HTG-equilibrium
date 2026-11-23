"""Compare and Queue tabs (WP7) with NiceGUI's simulated user (no Julia)."""
import asyncio
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from test_plot_core import make_tables
from test_results_tab import StubState

import htg
from htg import plotbuild as pb
from htg.gui.compare_tab import CAPTION, CompareTab
from htg.gui.differences import setting_differences, setting_label
from htg.gui.queue_tab import (
    EMPTY,
    WORKERS_CAPTION,
    QueueTab,
    point_text,
    run_text,
    worker_rows,
)
from htg.gui.runs import RunTables
from htg.runner import _current_point
from htg.settings import Settings

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = htg.load_case(ROOT / 'src/htg/templates/case_template.toml')


def run_tables(label, case, folder, comparison=None):
    parts = [make_tables(P=240.0), make_tables(P=300.0)]
    joined = {k: pd.concat([getattr(p, k) for p in parts], ignore_index=True) for k in ('summary', 'species', 'phases')}
    return RunTables(label, case, joined['summary'], joined['species'], joined['phases'], comparison, folder)


# ---------------------------------------------------------------- settings that differ
def test_setting_labels_follow_the_case_tab():
    other = TEMPLATE.replace(name='other', pressure_bar=300, solver='slsqp', multistart=3,
                             gas_carbon_min_fraction=0.9)
    table = setting_differences({'a': TEMPLATE, 'b': other, 'missing': None})
    assert list(table.columns) == ['Setting', 'a', 'b']
    rows = table.set_index('Setting')
    assert set(rows.index) == {'Pressures [bar]', 'Method', 'Initial guesses', 'Carbon split: gas, at least [%]'}
    assert rows.loc['Initial guesses', 'a'] == '(default)' and rows.loc['Initial guesses', 'b'] == '3'
    assert rows.loc['Carbon split: gas, at least [%]', 'b'] == '90'
    in_mpa = setting_differences({'a': TEMPLATE, 'b': other}, 'K', 'MPa').set_index('Setting')
    assert list(in_mpa.loc['Pressures [MPa]']) == ['24', '30']
    # version-2 keys
    v2 = {'a': {'solver': {'starts': 1}, 'feed': {'total_kg': 10.0, 'solids': {'value': 15.0, 'unit': '%'}},
                'constraints': {'carbon': {'gas': {'min': 0.9}}}, 'aqueous': {'cutoff_C': 390.0}},
          'b': {'solver': {'starts': 5}, 'feed': {'total_kg': 10.0, 'solids': {'value': 2.0, 'unit': 'kg'}},
                'constraints': {'carbon': {'gas': {'min': 0.8}}}, 'aqueous': {'cutoff_C': 380.0}}}
    rows = setting_differences(v2, 'K').set_index('Setting')
    assert list(rows.loc['Initial guesses']) == ['1', '5']
    assert list(rows.loc['Solids']) == ['15 %', '2 kg']
    assert list(rows.loc['Carbon split: gas, at least [%]']) == ['90', '80']
    assert list(rows.loc['Remove at and above [K]']) == ['663.15', '653.15']
    assert setting_label('phases.organic_liquid') == 'Possible products: Organic liquid'
    assert setting_label('constraints.minimum_amounts.gas.water') == 'Minimum amount: gas, water [mol / mol C]'
    assert setting_label('something.new') == 'something.new'


# ---------------------------------------------------------------- compare plot
def test_compare_spec(tmp_path):
    runs = [pb.CompareRun('a', run_tables('a', None, tmp_path), 0.2, 0.1),
            pb.CompareRun('b', run_tables('b', None, tmp_path))]
    items = pb.compare_items(runs)
    assert list(items)[0] == 'cge' and items['cge'] == 'Carbon gasification efficiency'
    assert items['gas:methane'] == 'Gas · methane'
    spec = pb.compare_runs(runs, 'gas:carbon dioxide', pb.PlotOptions(amount='g'))
    assert spec.title == 'Gas · carbon dioxide · g / kg dry feed'
    series = spec.panels[0].series
    assert [s.name for s in series] == ['a, 240 bar', 'a, 300 bar', 'b, 240 bar', 'b, 300 bar']
    assert [s.colour for s in series] == [0, 0, 1, 1] and [s.dash for s in series] == ['solid', 'dash'] * 2
    assert series[0].y[0] == pytest.approx(10.0 * 44.0)
    # total feed: runs without a total mass are left out
    spec = pb.compare_runs(runs, 'gas:carbon dioxide', pb.PlotOptions(basis='total'))
    assert [s.name for s in spec.panels[0].series] == ['a, 240 bar', 'a, 300 bar']
    assert spec.panels[0].series[0].y[0] == pytest.approx(10.0 * 0.2)
    cge = pb.compare_runs(runs, 'cge')
    assert cge.panels[0].y_label == 'Carbon gasification efficiency [%]' and len(cge.panels[0].series) == 4


def _simulate(page_factory, scenario):
    from nicegui.testing.user_simulation import user_simulation

    holder = {}

    def page():
        holder['tab'] = page_factory()

    async def main():
        async with user_simulation(page) as user:
            await user.open('/')
            await scenario(user, holder['tab'])

    asyncio.run(main())


def test_compare_tab(tmp_path, config_dir):
    measured = pd.DataFrame({'T_C': [300.0, 400.0], 'gas:carbon dioxide': [9.0, 21.0]})
    case_a = TEMPLATE.replace(name='a', comparison={'file': 'x.csv', 'quantity': 'mol_per_kg_dry_feed',
                                                    'label': 'Paper'})
    case_b = TEMPLATE.replace(name='b', multistart=4, pressure_bar=300)
    state = StubState({'saved:a': run_tables('a', case_a, tmp_path / 'a', measured),
                       'saved:b': run_tables('b', case_b, tmp_path / 'b')}, Settings.load(config_dir))

    async def scenario(user, tab):
        await user.should_see(CAPTION)
        await user.should_see('Choose two or more runs.')
        with user:
            tab.runs.set_value(['saved:a', 'saved:b'])
        await user.should_see('2 runs')
        assert tab.key.value == 'cge' and 'Gas · methane' in tab.key.options.values()
        assert tab.spec.title == 'Carbon gasification efficiency · %'
        settings = [row['Setting'] for row in tab.differences.rows]
        assert 'Initial guesses' in settings and 'Pressures [bar]' in settings and 'Measured data: label in the legend' in settings
        with user:
            tab.p_toggle.set_value('MPa')
            tab.t_toggle.set_value('K')
        settings = {row['Setting']: row for row in tab.differences.rows}
        assert settings['Pressures [MPa]']['b'] == '30'
        assert tab.spec.panels[0].x_label == 'Temperature [K]'
        assert tab.spec.panels[0].series[0].name == 'a, 24 MPa'
        assert tab.spec.panels[0].series[0].x[0] == pytest.approx(573.15)
        # measured data and the amount options
        assert tab.measured.options == {'': 'None', 'a': 'Paper (from a)'}
        assert not tab.amount.enabled
        with user:
            tab.key.set_value('gas:carbon dioxide')
            tab.measured.set_value('a')
        assert tab.amount.enabled
        marks = [s for s in tab.spec.panels[0].series if s.measured]
        assert marks and marks[0].name == 'Paper' and marks[0].y == [9.0, 21.0]
        assert tab.basis.props['options'][0]['disable'] is True  # the template cases have no total mass
        # Save plot…: one panel only
        user.find('Save plot…').click()
        dialog = tab.save_dialog
        assert list(dialog.layout.options) == ['one']
        with user:
            dialog.format.set_value('SVG')
        user.find('Save…').click()
        response = await user.download.next(timeout=10)
        assert b'<svg' in response.content[:2000]

    _simulate(lambda: CompareTab(state), scenario)


# ---------------------------------------------------------------- queue
def test_queue_wording():
    status = {'workers': [{'id': 2, 'state': 'busy', 'message': '', 'case': 'manure', 'point': (400.0, 240.0)},
                          {'id': 1, 'state': 'warming', 'message': 'compiling the equation of state', 'case': '',
                           'point': None}], 'queued_pieces': 0}
    rows = worker_rows(status)
    assert rows[0] == {'id': 1, 'state': 'warming', 'case': '', 'working_on': 'compiling the equation of state'}
    assert rows[1]['working_on'] == 'T = 400 °C, P = 240 bar'
    assert point_text((400.0, 240.0), 'K', 'MPa') == 'T = 673.15 K, P = 24 MPa'
    run = SimpleNamespace(done_points=8, total=13, failed_points=0, started=None, finished_running=False,
                          output_dir=None)
    assert run_text(run) == '8 of 13 points'
    task = SimpleNamespace(chain=[(300.0, 240.0), (350.0, 240.0)], run=SimpleNamespace(outcomes={(300.0, 240.0): 1}))
    assert _current_point(task) == (350.0, 240.0) and _current_point(None) is None


def test_queue_tab_page(tmp_path, config_dir):
    from nicegui import ui

    folder = tmp_path / 'results' / 'manure-1'
    folder.mkdir(parents=True)
    (folder / 'summary.csv').write_text('T_C\n')
    run = SimpleNamespace(id=3, name='manure', status='done', done_points=13, total=13, failed_points=1,
                          started=1.0, finished=241.0, finished_running=True, output_dir=folder, log=['started'])
    status = {'workers': [{'id': 1, 'state': 'busy', 'message': '', 'case': 'manure', 'point': (400.0, 240.0)}],
              'queued_pieces': 0}
    pool = SimpleNamespace(target=1, runs=[], status=lambda: status, cancel=lambda r: None)
    settings = Settings.load(config_dir)
    settings.set('pressure_unit', 'MPa')
    state = SimpleNamespace(pool=pool, settings=settings, results_dir=tmp_path / 'results')

    async def scenario(user, tab):
        await user.should_see(WORKERS_CAPTION)
        await user.should_see(EMPTY)
        assert [c['label'] for c in tab.worker_table.columns] == ['Worker', 'State', 'Case', 'Working on']
        assert tab.worker_table.rows[0]['working_on'] == 'T = 400 °C, P = 24 MPa'
        pool.runs.append(run)
        with user:
            tab.update()
        await user.should_see('#3 manure')
        await user.should_see('13 of 13 points, 1 failed · 4.0 min · saved in results/manure-1')
        await user.should_not_see(EMPTY)
        user.find(kind=ui.button, content='Results').click()

    shown = []
    _simulate(lambda: QueueTab(state, shown.append), scenario)
    assert shown == ['live:3']
