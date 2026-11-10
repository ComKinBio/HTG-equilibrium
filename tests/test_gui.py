"""The browser app: case-file editing, result tables and plots, and the page itself.

The page test uses NiceGUI's simulated user and stand-in workers (no Julia); the
plot tests use a small real T-P grid.
"""
import asyncio
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import htg
from htg.gui import plots
from htg.gui.app import State, add_species, build_page, table_rows
from htg.gui.document import CaseDocument, case_from_data
from htg.gui.editor import parse_elements, parse_numbers, parse_points
from htg.gui.runs import case_differences, load_saved, output_folder, saved_runs
from htg.runner import WorkerPool

ROOT = Path(__file__).resolve().parents[1]
CASES = sorted((ROOT / 'cases').glob('*.toml'))


# ---------------------------------------------------------------- text fields
def test_number_lists():
    assert parse_numbers('300, 350;400 450') == [300, 350, 400, 450]
    assert parse_numbers('300:400:25') == [300, 325, 350, 375, 400]
    assert parse_numbers('300:390:25') == [300, 325, 350, 375]
    assert parse_numbers('100, 300:350:50') == [100, 300, 350]
    assert parse_numbers('') == []
    with pytest.raises(ValueError):
        parse_numbers('300:400:0')
    with pytest.raises(ValueError):
        parse_numbers('hot')


def test_points_and_salt_elements():
    assert parse_points('300, 240\n\n350 300\n') == [[300, 240], [350, 300]]
    with pytest.raises(ValueError):
        parse_points('300')
    assert parse_elements('Na=1, H=1 C=1, O=3') == {'Na': 1, 'H': 1, 'C': 1, 'O': 3}
    assert parse_elements('Ca=0.5') == {'Ca': 0.5}
    with pytest.raises(ValueError):
        parse_elements('NaCl')


# ---------------------------------------------------------------- case documents
@pytest.mark.parametrize('path', CASES, ids=lambda p: p.stem)
def test_documents_match_case_files(path):
    """The editor reads a case file exactly like load_case, and an unchanged file saves unchanged."""
    document = CaseDocument.open(path)
    assert case_from_data(document.data, path) == htg.load_case(path)
    assert document.text() == path.read_text(encoding='utf-8')


def test_edits_keep_comments(tmp_path):
    path = tmp_path / 'case.toml'
    path.write_text(CaseDocument.new('edited').text(), encoding='utf-8')
    document = CaseDocument.open(path)
    comments = [line for line in path.read_text().splitlines() if line.strip().startswith('#')]
    document.data['feed']['water_g_per_g_dry'] = 7.5
    document.data['conditions'] = {'temperature_C': [300, 350], 'pressure_bar': [240, 300]}
    document.data.setdefault('constraints', {})['minimum_amounts'] = {'gas': {'methane': 0.01}}
    document.save()
    text = path.read_text()
    assert all(line in text for line in comments)
    case = htg.load_case(path)
    assert case.feed.water_g_per_g_dry == 7.5
    assert case.grid() == [(300, 240), (350, 240), (300, 300), (350, 300)]
    assert case.minimum_amounts == {('gas', 'methane'): 0.01}


def test_save_as_elsewhere_keeps_comparison_file(tmp_path):
    source = ROOT / 'cases' / 'yakaboylu2013_fig7_8.toml'
    document = CaseDocument.open(source)
    target = tmp_path / 'sub' / 'copy.toml'
    target.parent.mkdir()
    document.save(target)
    assert htg.load_case(target).comparison['file'] == htg.load_case(source).comparison['file']


def test_invalid_case_is_not_saved(tmp_path):
    document = CaseDocument.new('bad')
    document.data['conditions'] = {'temperature_C': 300, 'pressure_bar': -5}
    with pytest.raises(ValueError):
        document.save(tmp_path / 'bad.toml')
    assert not (tmp_path / 'bad.toml').exists()


def test_add_species_to_case_data():
    data = {'phases': {'gas': ['water']}}
    assert add_species(data, 'gas', 'methane') == 'the gas phase'
    assert data['phases']['gas'] == ['water', 'methane']
    with pytest.raises(ValueError):
        add_species(data, 'gas', 'methane')
    data['feed'] = {'elements_wt_percent': {'C': 50}}
    add_species(data, 'feed', 'cellulose')
    assert data['feed'] == {'species_moles': {'cellulose': 1.0}}


# ---------------------------------------------------------------- runs and tables
def test_case_differences():
    template = htg.load_case(ROOT / 'src/htg/templates/case_template.toml')
    other = template.replace(name='other', pressure_bar=300, solver='slsqp')
    table = case_differences({'a': template, 'b': other, 'missing': None})
    assert list(table.columns) == ['setting', 'a', 'b']
    assert set(table.setting) == {'conditions.pressure_bar', 'solver.method'}


def test_output_folders_do_not_collide(tmp_path):
    import datetime
    now = datetime.datetime(2026, 1, 2, 3, 4, 5)
    assert output_folder(tmp_path, 'x', now) == tmp_path / 'x-20260102-030405'


def test_table_rows_are_json_friendly():
    rows = table_rows(pd.DataFrame({'a': [1.234567891, np.nan], 'b': ['x', None]}))
    assert rows == [{'_row': 0, 'a': 1.23457, 'b': 'x'}, {'_row': 1, 'a': None, 'b': None}]


@pytest.fixture(scope='module')
def grid_run(base_model, tmp_path_factory):
    """A real 3 x 2 grid plus one failed point, saved and reloaded."""
    case = htg.load_case(ROOT / 'src/htg/templates/case_template.toml').replace(
        name='grid', temperature_C=[300, 350, 380], pressure_bar=[240, 300], aqueous_cutoff_C=390)
    points = []
    for chain in case.chains():
        htg.runner.solve_chain(base_model, case, chain, lambda T, P, outcome, start: points.append((T, P, outcome)))
    points.append((365.0, 240.0, 'made-up failure'))  # a failed point: a gap in the plots
    case = case.replace(temperature_C=[300, 350, 365, 380])
    sweep = htg.results.SweepResult(case, sorted(points, key=lambda p: (p[1], p[0])))
    folder = tmp_path_factory.mktemp('results') / 'grid-1'
    sweep.save(folder)
    return load_saved(folder, folder.parent)


def test_saved_runs_are_found(grid_run):
    assert saved_runs(grid_run.path.parent) == [grid_run.path]
    assert grid_run.pressures == [240, 300] and grid_run.temperatures == [300, 350, 365, 380]
    assert 'gas:methane' in grid_run.species_keys


def test_profiles_have_gaps_only_at_failed_points(grid_run):
    fig = plots.species_vs_temperature(grid_run.species, grid_run.summary, 'gas', log=False)
    assert len(fig.data) == 2 * len({t.legendgroup for t in fig.data})  # every species at both pressures
    for trace in fig.data:
        y = np.asarray(trace.y, dtype=float)
        if trace.name.endswith('240 bar'):
            assert list(trace.x) == [300, 350, 365, 380]
            assert np.isnan(y[2]) and not np.isnan(y[[0, 1, 3]]).any()
        else:  # 365 °C was only asked at 240 bar
            assert list(trace.x) == [300, 350, 380]
            assert not np.isnan(y).any()


def test_all_views_draw(grid_run):
    for group in plots.GROUPS:
        for quantity in plots.QUANTITIES:
            plots.species_vs_temperature(grid_run.species, grid_run.summary, group, quantity)
        plots.species_vs_pressure(grid_run.species, grid_run.summary, group, 350.0)
    for item in plots.map_items(grid_run):
        fig = plots.tp_map(grid_run, item)
        assert np.asarray(fig.data[0].z).shape == (2, 4)
    carbon = plots.carbon_vs_temperature(grid_run.phases, grid_run.summary)
    assert {t.legendgroup for t in carbon.data} >= {'gas'}
    both = plots.compare_runs([grid_run, grid_run], 'gas:methane')
    assert len(both.data) == 4


def test_gas_carbon_map_matches_summary(grid_run):
    fig = plots.tp_map(grid_run, 'gas_carbon_percent')
    z = np.asarray(fig.data[0].z, dtype=float)
    row = grid_run.summary[(grid_run.summary.T_C == 350) & (grid_run.summary.P_bar == 300)].iloc[0]
    assert z[1, 1] == pytest.approx(row.gas_carbon_percent)
    assert np.isnan(z[0, 2])  # the failed point


# ---------------------------------------------------------------- the page
def test_page_queues_a_case_and_shows_it(base_model, tmp_path):
    from nicegui import ui
    from nicegui.testing.user_simulation import user_simulation
    from test_runner import fake_worker

    async def scenario():
        with WorkerPool(2, worker_target=fake_worker) as pool:
            state = State(base_model, pool, results_dir=tmp_path / 'results', cases_dir=ROOT / 'cases')
            document = CaseDocument.new('page_test')
            document.data['conditions'] = {'temperature_C': [300, 350], 'pressure_bar': [240, 300]}
            state.document = document
            async with user_simulation(lambda: build_page(state)) as user:
                await user.open('/')
                await user.should_see('4 point(s) in 2 independent piece(s)')
                user.find('Add to queue').click()
                for _ in range(100):
                    if pool.runs and pool.runs[0].finished_running:
                        break
                    await asyncio.sleep(0.1)
                run = pool.runs[0]
                assert run.status == 'done' and run.total == 4
                await asyncio.sleep(1.5)  # the queue tab refreshes every second
                await user.should_see('4 of 4 points, 4 failed')
                user.find(kind=ui.button, content='Results').click()
                await user.should_see('0 of 4 points solved')
                assert (run.output_dir / 'summary.csv').exists()
                user.find('Search name or formula').type('methane')
                await user.should_see('Species')

    asyncio.run(scenario())
