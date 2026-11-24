"""WP4: the Case tab.

Pure functions (wording, notices, feed numbers) are tested directly; the form is
tested with NiceGUI's simulated user and a stand-in worker pool (no Julia).
"""
import asyncio
from pathlib import Path

import pytest

import htg
from htg.gui.app import State
from htg.gui.case_form import (check_feed, conflict_notice, element_error, entered_sum_notice, file_status, format_kg,
                               format_mol, initial_guesses_tooltip, parse_number, point_count, point_count_text,
                               problem_text)
from htg.gui.document import LEGACY_NOTICE, SALTS_NOTICE, CaseDocument, case_from_data
from htg.gui.editor import CaseEditor, parse_points
from htg.gui.species_tab import add_feed_species

ROOT = Path(__file__).resolve().parents[1]
FEED = {'total_kg': 10.0, 'solids': {'value': 15.0, 'unit': '%'}, 'ash': {'value': 8.45, 'unit': '%'},
        'elements_wt_percent': {'C': 30.4, 'H': 4.41, 'O': 53.58, 'N': 2.66, 'S': 0.51}}


class StubPool:
    def __init__(self):
        self.runs = []
        self.submitted = []

    def submit(self, case, output_dir):
        self.submitted.append(case)
        raise RuntimeError('not in this test')

    def status(self):
        return {'workers': [], 'queued_pieces': 0, 'target': 0}


# ---------------------------------------------------------------- pure functions
def test_feed_numbers():
    summary, field, message = check_feed(FEED)
    assert message is None
    assert [f'{summary["used"][e]:.2f}' for e in 'C H O N S'.split()] == ['33.20', '4.82', '58.52', '2.91', '0.56']
    rows = {row['name']: row for row in summary['rows']}
    assert (format_kg(rows['Water']['kg']), format_mol(rows['Water']['mol'])) == ('8.500', '471.8')
    assert format_kg(rows['Ash']['kg']) == '0.127'
    assert format_kg(rows['Ash-free solids']['kg']) == '1.373'
    assert (format_kg(rows['C']['kg']), format_mol(rows['C']['mol'])) == ('0.456', '37.96')
    assert entered_sum_notice(summary['entered_sum']) == (
        'The element percentages add up to 91.56 %. They have been scaled to 100 % (the Used column).')
    assert entered_sum_notice(100.0) is None and entered_sum_notice(99.999) is None


def test_feed_problems_name_their_field():
    feed = {**FEED, 'ash': {'value': 2, 'unit': 'kg'}}
    assert check_feed(feed)[1:] == ('ash', "Ash can't be more than the solids.")
    assert check_feed({**FEED, 'total_kg': 0})[1] == 'total'
    assert check_feed({**FEED, 'elements_wt_percent': {'C': -1}})[1] == 'element:C'
    assert check_feed({**FEED, 'elements_wt_percent': {}})[1] == 'elements'


def test_conflict_notice_conditions():
    data = {'conditions': {'temperature_C': [300, 400], 'pressure_bar': 240},
            'aqueous': {'cutoff_C': 390.0}, 'constraints': {'carbon': {'aqueous': {'exact': 0.03}}}}
    text = conflict_notice(data)
    assert text == ('Aqueous carbon is set to exactly 3 % (Carbon split). At 390 °C and above there is no '
                    'aqueous phase, so that setting is ignored there. Below 390 °C it applies.')
    assert 'At 663.15 K and above' in conflict_notice(data, 'K')
    data['constraints']['carbon']['aqueous'] = {'min': 0.05}
    assert conflict_notice(data).startswith('Aqueous carbon is set to at least 5 %')
    for change in ({'aqueous': {}},                                                     # no cut-off
                   {'constraints': {'carbon': {'aqueous': {'max': 0.05}}}},             # not exactly / at least
                   {'constraints': {'carbon': {'aqueous': {'exact': 0.0}}}},            # zero
                   {'conditions': {'temperature_C': [300, 380], 'pressure_bar': 240}}):  # all below
        assert conflict_notice({**data, **change}) is None


def test_initial_guesses_tooltip():
    assert initial_guesses_tooltip(1).startswith('You will start from 1 initial guess. The optimization')
    assert initial_guesses_tooltip(5).startswith('You will start from 5 initial guesses. The optimization')
    assert initial_guesses_tooltip(5).endswith('they are listed on the Results tab so you can compare them.')


def test_small_helpers():
    assert point_count({'temperature_C': [300, 350, 400], 'pressure_bar': [240, 300]}) == 6
    assert point_count({'temperature_C': [300, 400], 'pressure_bar': 240}, cutoff_C=390) == 3
    assert point_count({'points': [[300, 240], [350, 250]]}) == 2
    assert point_count({'temperature_C': [300]}) is None
    assert (point_count_text(13), point_count_text(1), point_count_text(None)) == ('13 points', '1 point', '')
    assert element_error('Cl') is None
    assert element_error('Qq') == '"Qq" is not an element symbol.'
    assert element_error('C', {'C': 1}) == 'C is already in the table.'
    assert parse_number(' 2.5 ') == 2.5 and parse_number('') is None
    with pytest.raises(ValueError, match='Give a value'):
        parse_number('', required=True)
    with pytest.raises(ValueError, match='^Initial guesses: give a whole number.$'):
        parse_number('1.5', 'Initial guesses', integer=True)
    assert problem_text('3 · Ash', "Ash can't be more than the solids.") == "Ash can't be more than the solids."
    assert problem_text('1 · Total mass [kg]', 'Give a value.') == 'Total mass: give a value.'
    assert parse_points('573.15, 24\n623.15 30', 'K', 'MPa') == [[300, 240], [350, 300]]
    assert file_status(None, True) == 'New case · not saved yet'
    assert file_status(ROOT / 'cases' / 'x.toml', True, ROOT) == 'cases/x.toml · unsaved changes'


def test_add_feed_species_writes_format_2():
    data = {'feed': {'total_kg': 10.0, 'solids': {'value': 15.0, 'unit': '%'},
                     'elements_wt_percent': {'C': 50}, 'salts': {'KCl': {'mol_per_g_dry_feed': 0.0,
                                                                         'elements': {'K': 1, 'Cl': 1}}}}}
    assert add_feed_species(data, 'cellulose') == 'the feed'
    assert data['feed']['species'] == {'water': {'value': 8.5, 'unit': 'kg'},
                                       'cellulose': {'value': 1.5, 'unit': 'kg'}}
    assert set(data['feed']) == {'species', 'salts'}
    add_feed_species(data, 'glucose')
    assert data['feed']['species']['glucose'] == {'value': 1.0, 'unit': 'kg'}
    with pytest.raises(ValueError, match='already in the feed'):
        add_feed_species(data, 'glucose')


# ---------------------------------------------------------------- the form
def run_form(model, scenario, document=None, cases_dir=None):
    """Run ``scenario(user, editor, state)`` on a page with the Case tab only."""
    from nicegui.testing.user_simulation import user_simulation
    state = State(model, StubPool(), results_dir='results', cases_dir=cases_dir or ROOT / 'cases')
    state.document = document or CaseDocument.new('form_test')
    holder = {}

    def page():
        holder['editor'] = CaseEditor(state, on_queued=lambda run: None)

    async def main():
        async with user_simulation(page) as user:
            await user.open('/')
            with user:
                await scenario(user, holder['editor'], state)

    asyncio.run(main())
    return state


def type_into(user, marker, text):
    user.find(marker=marker).clear().type(text)


def set_value(user, marker, value):
    for element in user.find(marker=marker).elements:
        element.value = value


def test_feed_card_numbers(base_model):
    async def scenario(user, editor, state):
        type_into(user, 'case-feed.total', '10')
        type_into(user, 'case-feed.solids', '15')
        type_into(user, 'case-feed.ash', '8.45')
        for element, value in zip('C H O N S'.split(), ('30.4', '4.41', '53.58', '2.66', '0.51')):
            type_into(user, f'case-feed.element:{element}', value)
        for element, used in zip('C H O N S'.split(), ('33.20', '4.82', '58.52', '2.91', '0.56')):
            assert editor.used_labels[element].text == used
        await user.should_see('91.56', marker='case-entered-sum')
        await user.should_see('The element percentages add up to 91.56 %. They have been scaled to 100 % '
                              '(the Used column).')
        rows = summary_rows(editor)
        assert rows['Water'] == ['8.500', '471.8']
        assert rows['Ash'][0] == '0.127' and rows['Ash-free solids'][0] == '1.373'
        assert rows['C'] == ['0.456', '37.96']
        feed = state.document.data['feed']
        assert feed['total_kg'] == 10 and feed['solids'] == {'value': 15, 'unit': '%'}
        assert feed['ash'] == {'value': 8.45, 'unit': '%'}
        assert editor.queue_button.enabled

        # ash as kg, more than the 1.5 kg of solids
        set_value(user, 'case-feed.ash-unit', 'kg')
        assert feed['ash']['unit'] == 'kg' and feed['ash']['value'] == pytest.approx(0.12675)
        type_into(user, 'case-feed.ash', '2')
        ash = user.find(marker='case-feed.ash').elements.pop()
        assert ash.error == "Ash can't be more than the solids."
        assert not editor.queue_button.enabled and not editor.save_button.enabled
        await user.should_see("✗ Ash can't be more than the solids.")
        type_into(user, 'case-feed.ash', '0.1')
        assert ash.error is None and editor.queue_button.enabled

        # adding and removing elements
        type_into(user, 'case-new-element', 'qq')
        user.find(marker='case-add-element').click()
        await user.should_see('"Qq" is not an element symbol.')
        type_into(user, 'case-new-element', 'cl')
        user.find(marker='case-add-element').click()
        assert list(feed['elements_wt_percent']) == ['C', 'H', 'O', 'N', 'S', 'Cl']
        user.find(marker='case-remove-S').click()
        assert 'S' not in feed['elements_wt_percent']

        # a field that is not a number blocks saving and names the field
        type_into(user, 'case-feed.total', 'ten')
        await user.should_see('✗ Total mass: "ten" is not a number.')
        assert feed['total_kg'] == 10 and not editor.save_button.enabled

    run_form(base_model, scenario)


def summary_rows(editor):
    """{row name: [kg, mol]} as shown in the Feed summary."""
    from nicegui import ui
    grid = next(iter(editor_elements(editor, 'case-feed-summary')))
    labels = [child.text for child in grid.default_slot.children if isinstance(child, ui.label)][3:]
    return {labels[i]: labels[i + 1:i + 3] for i in range(0, len(labels), 3)}


def editor_elements(editor, marker):
    from nicegui import ElementFilter
    return list(ElementFilter(marker=marker))


def test_conflict_notice_in_the_form(base_model):
    async def scenario(user, editor, state):
        data = state.document.data
        type_into(user, 'case-conditions.temperatures', '300, 400')
        type_into(user, 'case-aqueous.cutoff', '390')
        set_value(user, 'case-carbon-aqueous', 'exact')
        assert not editor.conflict_label.visible
        type_into(user, 'case-carbon.aqueous', '3')
        assert data['constraints']['carbon']['aqueous'] == {'exact': 0.03}
        await user.should_see('Aqueous carbon is set to exactly 3 % (Carbon split). At 390 °C and above')
        assert editor.conflict_label.visible
        # 1: the cut-off
        type_into(user, 'case-aqueous.cutoff', '')
        assert not editor.conflict_label.visible
        type_into(user, 'case-aqueous.cutoff', '390')
        assert editor.conflict_label.visible
        # 2: the aqueous carbon choice
        set_value(user, 'case-carbon-aqueous', 'max')
        assert not editor.conflict_label.visible and data['constraints']['carbon']['aqueous'] == {'max': 0.03}
        set_value(user, 'case-carbon-aqueous', 'min')
        assert editor.conflict_label.visible
        await user.should_see('Aqueous carbon is set to at least 3 %')
        type_into(user, 'case-carbon.aqueous', '0')
        assert not editor.conflict_label.visible
        type_into(user, 'case-carbon.aqueous', '3')
        assert editor.conflict_label.visible
        # 3: a temperature at or above the cut-off
        type_into(user, 'case-conditions.temperatures', '300, 380')
        assert not editor.conflict_label.visible
        type_into(user, 'case-conditions.temperatures', '300:400:50')
        assert editor.conflict_label.visible
        # the notice and the fields follow the temperature unit
        set_value(user, 'case-temperature-unit', 'K')
        assert 'At 663.15 K and above' in editor.conflict_label.text
        assert user.find(marker='case-conditions.temperatures').elements.pop().value == '573.15, 623.15, 673.15'
        assert state.settings.get('temperature_unit') == 'K'
        await user.should_see('Temperatures [K]')
        await user.should_see('Remove at and above [K]')
        type_into(user, 'case-conditions.temperatures', '573.15, 773.15')
        assert data['conditions']['temperature_C'] == [300, 500]

    run_form(base_model, scenario)


def test_initial_guesses_tooltip_in_the_form(base_model):
    async def scenario(user, editor, state):
        type_into(user, 'case-solver.starts', '1')
        assert editor.starts_tooltip.text.startswith('You will start from 1 initial guess.')
        type_into(user, 'case-solver.starts', '7')
        assert editor.starts_tooltip.text.startswith('You will start from 7 initial guesses.')
        assert state.document.data['solver']['starts'] == 7
        type_into(user, 'case-solver.starts', '0')
        await user.should_see('✗ Initial guesses: give a value of 1 or more.')
        assert state.document.data['solver']['starts'] == 7

    run_form(base_model, scenario)


def test_conditions_and_point_count(base_model, config_dir):
    async def scenario(user, editor, state):
        data = state.document.data
        type_into(user, 'case-conditions.temperatures', '300:600 n=13')
        await user.should_see('14 points')  # the template's cut-off (390 °C) adds a point
        type_into(user, 'case-aqueous.cutoff', '')
        await user.should_see('13 points')
        set_value(user, 'case-pressure-unit', 'MPa')
        type_into(user, 'case-conditions.pressures', '24, 30')
        assert data['conditions'] == {'temperature_C': [300 + 25 * i for i in range(13)], 'pressure_bar': [240, 300]}
        await user.should_see('26 points')
        set_value(user, 'case-conditions-mode', 'points')
        assert len(data['conditions']['points']) == 26 and 'temperature_C' not in data['conditions']
        type_into(user, 'case-conditions.points', '300, 24\n350, 25')
        assert data['conditions'] == {'points': [[300, 240], [350, 250]]}
        await user.should_see('2 points')
        user.find(marker='case-input-formats').elements.pop().value = True
        assert state.settings.get('input_formats_open') is True

    run_form(base_model, scenario)
    settings = htg.settings.Settings.load(config_dir)
    assert settings.get('pressure_unit') == 'MPa' and settings.get('input_formats_open') is True


def test_species_amounts(base_model):
    async def scenario(user, editor, state):
        set_value(user, 'case-feed-mode', 'species')
        feed = state.document.data['feed']
        assert feed == {'species': {'water': {'value': 4.5, 'unit': 'kg'}, 'cellulose': {'value': 0.5, 'unit': 'kg'}}}
        await user.should_see('Water (always present)')
        await user.should_see('Solids (not water)')
        type_into(user, 'case-feed.water', '9')
        assert feed['species']['water'] == {'value': 9, 'unit': 'kg'}
        set_value(user, 'case-feed.water-unit', 'mol')
        assert feed['species']['water']['unit'] == 'mol'
        assert feed['species']['water']['value'] == pytest.approx(9000 / 18.015, rel=1e-5)
        user.find(marker='case-add-species').click()
        await user.should_see('✗ Species: choose a species or remove the row.')
        set_value(user, 'case-species-1', 'cellulose')
        await user.should_see('cellulose is already in the feed.')
        set_value(user, 'case-species-1', 'methanol')
        assert list(state.document.data['feed']['species']) == ['water', 'cellulose', 'methanol']
        assert editor.queue_button.enabled
        case = case_from_data(state.document.data)
        assert case.feed_v2['species']['methanol'] == {'value': 1.0, 'unit': 'kg'}
        set_value(user, 'case-feed-mode', 'elements')
        assert 'total_kg' in state.document.data['feed']
        set_value(user, 'case-feed-mode', 'species')  # the species come back
        assert 'methanol' in state.document.data['feed']['species']

    run_form(base_model, scenario)


def test_save_writes_format_2(base_model, tmp_path):
    path = tmp_path / 'old.toml'
    path.write_text((ROOT / 'tests/data/legacy_cases/manure_gas_min.toml').read_text(encoding='utf-8'))

    async def scenario(user, editor, state):
        editor.open(path)
        await user.should_see(LEGACY_NOTICE)
        type_into(user, 'case-solver.starts', '3')
        await user.should_see('old.toml · unsaved changes')
        editor._save()
        await user.should_see('old.toml', marker=None)
        assert not editor.dirty()
        assert not any(LEGACY_NOTICE in e.text for e in editor_elements(editor, 'case-notice'))

    run_form(base_model, scenario)
    text = path.read_text()
    assert 'format = 2' in text and 'starts = 3' in text and 'multistart' not in text
    assert '[constraints.carbon]' in text or 'carbon =' in text
    assert htg.load_case(path).starts == 3


def test_salts_and_blocked_files(base_model, tmp_path):
    async def salts(user, editor, state):
        editor.open(ROOT / 'cases' / 'manure_ions_minerals.toml')
        await user.should_see(SALTS_NOTICE)
        assert editor.queue_button.enabled

    run_form(base_model, salts)

    blocked = tmp_path / 'unmodelled.toml'
    blocked.write_text('name = "old"\n[conditions]\ntemperature_C = 300\npressure_bar = 240\n'
                       '[feed]\nelements_wt_percent = { C = 50, H = 6, O = 40 }\nunmodelled_wt_percent = 4\n'
                       'water_g_per_g_dry = 4\n[phases]\ngas = ["water", "methane", "carbon dioxide"]\n',
                       encoding='utf-8')

    async def unmodelled(user, editor, state):
        editor.open(blocked)
        await user.should_see("This case uses unmodelled_wt_percent, which the form can't show. "
                              'Edit the file directly.')
        assert not editor.fields and editor.queue_button.enabled
        await user.should_see('Ready to run: 1 point in 1 independent piece.')

    run_form(base_model, unmodelled)


def test_check_button(base_model):
    async def scenario(user, editor, state):
        user.find('Check').click()
        await user.should_see('All 5 points passed the checks.', retries=200)
        set_value(user, 'case-carbon-gas', 'min')
        type_into(user, 'case-carbon.gas', '90')
        set_value(user, 'case-carbon-solid', 'min')
        type_into(user, 'case-carbon.solid', '20')
        user.find('Check').click()
        await user.should_see('The carbon shares add up to more than 100 %.', retries=200)

    run_form(base_model, scenario)
