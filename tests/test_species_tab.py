"""WP8: the Species tab, the new / edit species dialog, import, and restarting the workers.

The page tests use NiceGUI's simulated user, a stand-in worker pool and a user
species database in tmp_path (never the real user folder).
"""
import asyncio
import time

import pandas as pd
import pytest

import htg
from htg.gui.app import State
from htg.gui.editor import CaseEditor
from htg.gui.species_dialog import format_number, parse_number
from htg.gui.species_tab import SpeciesTab, add_species, element_error
from htg.runner import WorkerPool
from htg.species_catalogue import count_text, filter_catalogue
from htg.user_species import UserDatabase, format_nasa_record, spec_from_model
from test_user_species import DATA, TEMPLATE

T = 623.15


@pytest.fixture(autouse=True)
def no_real_user_folder(tmp_path, monkeypatch):
    monkeypatch.setenv('HTG_USER_DATA', str(tmp_path / 'env-user-folder'))


class StubPool:
    """Stands in for WorkerPool: records restarts."""

    def __init__(self, active=False):
        self.runs = []
        self.restarts = []
        self.active = active

    def restart(self, model_args):
        self.restarts.append(dict(model_args))
        return not self.active

    def status(self):
        return {'workers': [], 'queued_pieces': 0, 'target': 0}


def methane_m(model):
    return model.registry.resolve_parameters(model.databases, T, ['methane'])['methane']['m']


def run_page(model, scenario, pool=None):
    """Run ``scenario(user, tab, state)`` on a page with the case editor and the Species tab."""
    from nicegui.testing.user_simulation import user_simulation
    state = State(model, pool or StubPool(), results_dir='results', cases_dir='cases')
    holder = {}

    def page():
        editor = CaseEditor(state, on_queued=lambda run: None)
        holder['tab'] = SpeciesTab(state, editor)

    async def main():
        async with user_simulation(page) as user:
            await user.open('/')
            with user:  # direct calls on the tab create UI too
                await scenario(user, holder['tab'], state)

    asyncio.run(main())
    return state


async def settle(condition, timeout=10.0):
    """Wait until condition() is true (async handlers run as background tasks)."""
    end = time.time() + timeout
    while not condition():
        if time.time() > end:
            raise AssertionError('condition not reached')
        await asyncio.sleep(0.05)


# ---------------------------------------------------------------- small helpers
def test_number_fields():
    assert parse_number(' 1.5E-3 ') == 1.5e-3 and parse_number('2.0D+01') == 20.0
    assert parse_number('') is None and parse_number('abc') is None and parse_number('nan') is None
    assert format_number(None) == '' and format_number(300.0) == '300'
    value = 0.12345678901234567
    assert parse_number(format_number(value)) == value
    assert element_error('C, H, O') is None
    assert element_error('C, Qq') == 'Unknown element: Qq. Use symbols such as C, H, O.'


def test_add_species_to_phases_and_feed():
    data = {'phases': {'gas': ['water']}}
    assert add_species(data, 'aqueous', 'Na+') == 'the aqueous phase'
    assert data['phases']['aqueous'] == ['Na+']
    data['feed'] = {'elements_wt_percent': {'C': 50}}
    assert add_species(data, 'feed', 'cellulose') == 'the feed'
    assert data['feed'] == {'species': {'water': {'value': 0.0, 'unit': 'kg'},
                                        'cellulose': {'value': 1.0, 'unit': 'kg'}}}


# ---------------------------------------------------------------- restarting the workers
def tagging_worker(worker_id, inbox, outbox, model_args, warm_up):
    """A stand-in worker whose outcomes name the model arguments it was started with."""
    outbox.put(('state', worker_id, 'idle', ''))
    while True:
        message = inbox.get()
        if message is None:
            break
        task_id, case, chain = message
        for T_, P in chain:
            time.sleep(0.3)
            outbox.put(('point', worker_id, task_id, T_, P, f'data {model_args.get("tag")}', 'test'))
        outbox.put(('task_done', worker_id, task_id))
    outbox.put(('state', worker_id, 'stopped', ''))


def test_pool_restart_waits_for_submitted_runs(tmp_path):
    case = htg.load_case(TEMPLATE).replace(name='t', temperature_C=[300, 350], pressure_bar=240,
                                           continuation=True)
    with WorkerPool(1, {'tag': 'old'}, worker_target=tagging_worker) as pool:
        first = pool.submit(case)
        time.sleep(0.2)
        assert pool.restart({'tag': 'new'}) is False  # a run is going: it keeps the old data
        second = pool.submit(case.replace(name='later'))
        assert first.wait(20) and second.wait(20)
        assert {outcome for outcome, _ in first.outcomes.values()} == {'data old'}
        assert {outcome for outcome, _ in second.outcomes.values()} == {'data new'}
        assert pool.model_args == {'tag': 'new'}
        # Nothing running: the workers are replaced at once.
        assert pool.restart({'tag': 'newer'}) is True
        third = pool.submit(case.replace(name='third'))
        assert third.wait(20)
        assert {outcome for outcome, _ in third.outcomes.values()} == {'data newer'}
        for _ in range(50):
            if sum(w['state'] != 'stopped' for w in pool.status()['workers']) == 1:
                break
            time.sleep(0.1)
        assert sum(w['state'] != 'stopped' for w in pool.status()['workers']) == 1


# ---------------------------------------------------------------- the tab
def test_filters_count_and_add_buttons(tmp_path):
    model = htg.Model.load(user_data=tmp_path / 'user')

    async def scenario(user, tab, state):
        assert tab.count.text == count_text(filter_catalogue(tab.catalogue), tab.catalogue)
        user.find('Contains elements').type('C, H')
        user.find('Exact match').click()
        expected = filter_catalogue(tab.catalogue, elements='C, H', exact=True)
        assert tab.count.text == count_text(expected, tab.catalogue)
        assert {row['Species'] for row in tab.table.rows} == set(expected.Species)
        tab.phase_chips['solid'].set_selected(True)
        solids = filter_catalogue(tab.catalogue, elements='C, H', exact=True, can_be_in=['solid'])
        assert tab.count.text == count_text(solids, tab.catalogue) and len(solids) < len(expected)
        tab.phase_chips['solid'].set_selected(False)
        user.find('Contains elements').type(', Qq')  # an unknown symbol is ignored
        assert tab.count.text == count_text(expected, tab.catalogue)

        tab.select('methane')
        assert tab.selected_label.text == 'methane (CH4)'
        enabled = {t for t, b in tab.add_buttons.items() if b.enabled}
        assert enabled == {'gas', 'aqueous', 'organic_liquid', 'feed'}
        tab.select('Na+')
        assert {t for t, b in tab.add_buttons.items() if b.enabled} == {'aqueous'}
        tab.add('aqueous')
        assert 'Na+' in state.document.data['phases']['aqueous']
        await user.should_see('Adds it to the case open on the Case tab')
        await user.should_see('New and edited species are saved in your own database; '
                              'the shipped databases are never changed.')
        assert not tab.remove_button.visible  # not a user species

    run_page(model, scenario)


def test_database_off_disables_editing():
    model = htg.Model.load(user_data=None)

    async def scenario(user, tab, state):
        await user.should_see('Your own species database is off (--no-user-species).')
        assert not tab.new_button.enabled and not tab.import_button.enabled
        tab.select('methane')
        assert not tab.edit_button.enabled

    run_page(model, scenario)


def test_new_species_dialog_saves_and_reloads(tmp_path):
    folder = tmp_path / 'user'
    model = htg.Model.load(user_data=folder)
    r = model.databases.nasa['CH4']
    record = format_nasa_record('TESTANOL', {'C': 9, 'H': 12, 'O': 3}, 'G', r['T_low'], r['T_mid'],
                                r['T_high'], r['high'], r['low'])
    pool = StubPool()

    async def scenario(user, tab, state):
        user.find('+ New species').click()
        dialog = tab.dialog
        await user.should_see('Thermochemistry')
        user.find(marker='species-name').type('Testanol')
        user.find(marker='species-formula').type('C9H12O3')
        await user.should_see('C 9 · H 12 · O 3 · 168.19 g/mol')
        await user.should_see('Formula understood')
        dialog.thermo_radio.set_value('paste')
        user.find(marker='nasa-paste').type(record)
        await user.should_see('Thermochemistry complete')
        await user.should_see('Missing equation of state parameters: m, σ, ε/k')
        await user.should_see('Can be in: gas, aqueous, organic liquid')
        user.find(marker='eos-m').type('3.1')
        user.find(marker='eos-sigma').type('3.6')
        user.find(marker='eos-epsilon').type('250')
        await user.should_see('Equation of state parameters complete')
        assert dialog.check().usable
        user.find('Save species').click()
        await settle(lambda: 'testanol' in set(tab.catalogue.Species))
        row = tab.catalogue.set_index('Species').loc['testanol']
        assert row.Usable and row.Yours
        assert state.model.user_info['testanol']
        assert 'testanol' in state.phase_options['gas'] and 'testanol' in state.feed_species
        assert pool.restarts == [state.model.load_args]
        assert tab.choice['Species'] == 'testanol' and tab.remove_button.visible

        # Remove it again.
        user.find('Remove from my database').click()
        user.find(marker='confirm-remove').click()
        await settle(lambda: 'testanol' not in set(tab.catalogue.Species))
        assert len(pool.restarts) == 2 and tab.choice is None

    run_page(model, scenario, pool)
    assert 'testanol' not in UserDatabase.load(folder).owners()


def test_incomplete_species_is_saved(tmp_path):
    folder = tmp_path / 'user'
    model = htg.Model.load(user_data=folder)

    async def scenario(user, tab, state):
        tab.new_species()
        user.find(marker='species-name').type('halfdone')
        user.find(marker='species-formula').type('C2H6O2')
        assert not tab.dialog.check().usable
        user.find('Save species').click()
        await settle(lambda: 'halfdone' in set(tab.catalogue.Species))
        row = tab.catalogue.set_index('Species').loc['halfdone']
        assert not row.Usable and row.Yours

    run_page(model, scenario)


def test_edit_shipped_species_and_reset(tmp_path):
    folder = tmp_path / 'user'
    model = htg.Model.load(user_data=folder)
    original = methane_m(model)
    pool = StubPool(active=True)

    async def scenario(user, tab, state):
        tab.select('methane')
        user.find('Edit selected').click()
        await user.should_see('Edit methane')
        await user.should_not_see('Reset to original')
        assert tab.dialog.spec.m == original
        field = user.find(marker='eos-m')
        field.clear().type(format_number(original * 1.5))
        assert tab.dialog.check().usable
        user.find('Save species').click()
        await settle(lambda: tab.catalogue.set_index('Species').loc['methane'].Overrides)
        assert methane_m(state.model) == pytest.approx(original * 1.5)
        await user.should_see('Runs already queued use the old species data. Later runs use the new data.')

        user.find('Edit selected').click()
        await user.should_see('Reset to original')
        user.find('Reset to original').click()
        await settle(lambda: not tab.catalogue.set_index('Species').loc['methane'].Overrides)
        assert methane_m(state.model) == original

    run_page(model, scenario, pool)


def test_ion_dialog_charge_and_records(tmp_path):
    model = htg.Model.load(user_data=tmp_path / 'user')

    async def scenario(user, tab, state):
        tab.new_species()
        dialog = tab.dialog
        dialog.kind_toggle.set_value('ion')
        assert dialog.mode == 'hkf_record' and dialog.spec.model == 'ePC-SAFT'
        assert not dialog.model_select.enabled and not dialog.assoc_box.visible
        dialog.charge.set_value(1)
        assert 'Ag+' in dialog.hkf_select.options
        dialog.hkf_select.set_value('Ag+')
        await user.should_see('✓ Charge matches')
        dialog.kind_toggle.set_value('solid')
        assert dialog.mode == 'nasa' and not dialog.eos_box.visible
        dialog.thermo_radio.set_value('maier_kelley')
        await user.should_see('Single temperature range only. Minerals with phase transitions are not '
                              'supported yet.')

    run_page(model, scenario)


def _pcsaft_csv():
    ship = pd.read_csv(DATA / 'PROJECT_PCSAFT.csv')
    base = ship[(ship.Species == 'methane') & (ship.Model == 'PC-SAFT')].head(1)
    new = base.assign(Species='importane', Formula='C13H28O', m=4.2)
    existing = base.assign(m=1.234)
    broken = base.assign(Species='brokane', m='not a number')
    return pd.concat([new, existing, broken]).to_csv(index=False)


def test_import_dialog(tmp_path):
    folder = tmp_path / 'user'
    model = htg.Model.load(user_data=folder)
    pool = StubPool()

    async def scenario(user, tab, state):
        user.find('Import from file…').click()
        dialog = tab.dialog
        await user.should_see('Species that replace existing values are unticked by default. Imported '
                              'species go into your own database and can be removed again on the Species tab.')
        await dialog.load_text(_pcsaft_csv(), 'mine.csv')
        assert dialog.format.value == 'pcsaft'
        assert [row['What happens'] for row in dialog.table.rows][:2] == ['new', 'replaces existing values']
        assert dialog.table.rows[2]['What happens'].startswith('skipped: ')
        assert [row['Species'] for row in dialog.table.selected] == ['importane']
        # Skipped rows cannot be ticked.
        dialog.table.selected = list(dialog.table.rows)
        dialog.selection_changed()
        assert len(dialog.table.selected) == 2
        dialog.table.selected = [dialog.table.rows[0]]
        dialog.selection_changed()
        assert dialog.import_button.text == 'Import 1 species'
        user.find('Import 1 species').click()
        await settle(lambda: 'importane' in set(tab.catalogue.Species))
        assert tab.catalogue.set_index('Species').loc['importane'].Yours
        assert methane_m(state.model) == methane_m(model)  # the replacement was not ticked
        assert len(pool.restarts) == 1

    run_page(model, scenario, pool)


def test_unrecognised_import_file(tmp_path):
    model = htg.Model.load(user_data=tmp_path / 'user')

    async def scenario(user, tab, state):
        tab.import_file()
        await tab.dialog.load_text('just some text\n', 'notes.txt')
        assert tab.dialog.format.value is None
        await user.should_see('Format not recognised. Choose the format.')
        assert not tab.dialog.import_button.enabled

    run_page(model, scenario)


def test_edit_dialog_fills_existing_species(tmp_path):
    model = htg.Model.load(user_data=tmp_path / 'user')
    spec = spec_from_model(model, 'methane')
    assert spec.thermo == 'nasa_record'

    async def scenario(user, tab, state):
        tab.select('methane')
        tab.edit_selected()
        dialog = tab.dialog
        assert dialog.mode == 'nasa_record' and dialog.nasa_select.value == spec.nasa_record
        await user.should_see('Only record with formula CH4')
        assert not dialog.kind_toggle.enabled

    run_page(model, scenario)
