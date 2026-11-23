"""Case files in format 2, the feed arithmetic, and the conversion of old files.

The old (format 1) versions of the example cases are kept in
tests/data/legacy_cases; old and converted forms must give the same Case in
every number the model uses.
"""
import copy
import math
import os
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import htg
from htg.case import (LegacyConversionError, carbon_bound_messages, case_from_dict, check_carbon_bounds,
                      effective_carbon_bounds, fluid_phase_count, is_v2, legacy_carbon_bounds, legacy_to_v2,
                      starts_from_multistart)
from htg.feed_v2 import (FeedError, UNMODELLED_MESSAGE, feed_fractions, feed_summary, feed_v2_to_model,
                         legacy_feed_to_v2, molar_mass_function)
from htg.gui.document import LEGACY_NOTICE, SALTS_NOTICE, CaseDocument, clean

ROOT = Path(__file__).resolve().parents[1]
CASES = ROOT / 'cases'
LEGACY = Path(__file__).resolve().parent / 'data' / 'legacy_cases'
TEMPLATE = ROOT / 'src' / 'htg' / 'templates' / 'case_template.toml'
LEGACY_FILES = sorted(LEGACY.glob('*.toml'))
# Elements of a few species, for species-amount feeds without the databases.
ELEMENTS = {'water': {'H': 2, 'O': 1}, 'cellulose': {'C': 6, 'H': 10, 'O': 5},
            'methanol': {'C': 1, 'H': 4, 'O': 1}}


def species_elements(name):
    return ELEMENTS[name]


def new_form(legacy_path):
    """The shipped format-2 file for an old file (the template lives elsewhere)."""
    return TEMPLATE if legacy_path.stem == 'case_template' else CASES / legacy_path.name


def model_numbers(case):
    """Everything the calculation takes from the case, as plain values."""
    feed = case.feed
    dry = feed.dry_element_wt_percent(None)
    return {'dry': dry, 'water': feed.total_water_g_per_g_dry(), 'ash_rest': 100 - sum(dry.values()),
            'salts': feed.salts, 'temperatures': case.temperatures(), 'grid': case.grid(),
            'chains': case.chains(), 'carbon_bounds': case.carbon_bounds,
            'solid_carbon_target': case.solid_carbon_target, 'starts': case.starts}


def same_numbers(a, b):
    assert a['dry'].keys() == b['dry'].keys()
    for element in a['dry']:
        assert a['dry'][element] == pytest.approx(b['dry'][element], rel=1e-12, abs=1e-12), element
    assert a['water'] == pytest.approx(b['water'], rel=1e-12)
    assert a['ash_rest'] == pytest.approx(b['ash_rest'], rel=1e-12, abs=1e-12)
    for key in ('salts', 'temperatures', 'grid', 'chains', 'carbon_bounds', 'solid_carbon_target', 'starts'):
        assert a[key] == b[key], key


def rest_of_case(case):
    """The case without the feed and temperature list (compared by number elsewhere)."""
    return case.replace(feed=htg.Feed(), feed_v2=None, temperature_C=None, refine=None, comparison=None,
                        description='')


# ---------------------------------------------------------------- old and converted files
@pytest.mark.parametrize('path', LEGACY_FILES, ids=lambda p: p.stem)
def test_old_and_converted_files_give_the_same_case(path):
    old = htg.load_case(path)
    shipped = htg.load_case(new_form(path))
    converted = case_from_dict(legacy_to_v2(tomllib.loads(path.read_text())), name=path.stem)
    for new in (shipped, converted):
        assert new.feed_v2 is not None
        same_numbers(model_numbers(old), model_numbers(new))
        assert rest_of_case(new) == rest_of_case(old)
    if old.feed.ash_wt_percent is not None:  # ash given in the old file is kept
        assert converted.feed.ash_wt_percent == old.feed.ash_wt_percent


@pytest.mark.parametrize('path', sorted(CASES.glob('*.toml')) + [TEMPLATE], ids=lambda p: p.stem)
def test_v2_round_trip(path, tmp_path):
    case = htg.load_case(path)
    assert is_v2(tomllib.loads(path.read_text()))
    text = case.to_toml()
    assert text.startswith('format = 2\n')
    (tmp_path / path.name).write_text(text)
    again = htg.load_case(tmp_path / path.name)
    again.comparison = case.comparison
    assert again == case
    assert again.to_toml() == text


def test_old_python_case_is_written_in_format_2(tmp_path):
    old = htg.load_case(LEGACY / 'manure_configurable.toml')
    old.save(tmp_path / 'c.toml')
    data = tomllib.loads((tmp_path / 'c.toml').read_text())
    assert data['format'] == 2 and 'refine' not in data.get('sweep', {})
    assert set(data['feed']) == {'total_kg', 'solids', 'ash', 'elements_wt_percent', 'salts'}
    same_numbers(model_numbers(old), model_numbers(htg.load_case(tmp_path / 'c.toml')))


def test_documents_convert_old_files_and_keep_comments(tmp_path):
    path = tmp_path / 'old.toml'
    path.write_text((LEGACY / 'manure_ions_minerals.toml').read_text())
    document = CaseDocument.open(path)
    assert document.converted_from_legacy and document.legacy_blocker is None
    assert document.has_salts and document.notices() == [LEGACY_NOTICE, SALTS_NOTICE]
    same_numbers(model_numbers(htg.load_case(path)), model_numbers(document.case()))
    comments = [line for line in path.read_text().splitlines()
                if line.startswith('#') and 'water_g_per_g_dry' not in line]
    document.save()
    assert not document.converted_from_legacy
    text = path.read_text()
    assert all(line in text for line in comments)
    data = tomllib.loads(text)
    assert data['format'] == 2
    assert 'water_g_per_g_dry' not in data['feed'] and 'gas_carbon_fraction' not in data['constraints']
    assert data['constraints']['carbon'] == {'gas': {'exact': 0.8}}
    assert set(data['feed']['salts']) == {'NaHCO3', 'KCl', 'CsCl', 'LiCl', 'MgCl2'}
    reopened = CaseDocument.open(path)
    assert not reopened.converted_from_legacy
    assert reopened.case() == document.case()


def test_document_keeps_files_that_cannot_be_converted(tmp_path):
    path = tmp_path / 'odd.toml'
    text = (LEGACY / 'case_template.toml').read_text().replace(
        'water_g_per_g_dry = 4.0', 'water_g_per_g_dry = 4.0\nunmodelled_wt_percent = 3.0\nash_wt_percent = 5.44')
    path.write_text(text)
    document = CaseDocument.open(path)
    assert not document.converted_from_legacy
    assert document.legacy_blocker == UNMODELLED_MESSAGE
    assert document.notices() == [UNMODELLED_MESSAGE]
    document.case()  # still loads
    document.save()
    assert tomllib.loads(path.read_text()) == clean(tomllib.loads(text))  # still the old format


# ---------------------------------------------------------------- legacy_to_v2 rules
def legacy(**changes):
    data = tomllib.loads((LEGACY / 'case_template.toml').read_text())
    for section, values in changes.items():
        if values is None:
            data.pop(section, None)
        else:
            data[section] = values
    return data


def test_legacy_feed_rules():
    dry = legacy_feed_to_v2({'elements_wt_percent': {'C': 45.0, 'H': 5.0, 'O': 40.0}, 'ash_wt_percent': 10.0,
                             'water_g_per_g_dry': 3.0, 'moisture_wt_percent': 20.0})
    assert dry['total_kg'] == pytest.approx(1 + 3.0 + 20 / 80)
    assert dry['solids'] == {'value': 1.0, 'unit': 'kg'}
    assert dry['ash'] == {'value': 10.0, 'unit': '%'}
    assert dry['elements_wt_percent'] == pytest.approx({'C': 50.0, 'H': 50 / 9, 'O': 400 / 9})
    daf = legacy_feed_to_v2({'elements_wt_percent': {'C': 50.0, 'H': 6.0, 'O': 44.0}, 'basis': 'dry_ash_free',
                             'ash_wt_percent': 5.0})
    assert daf['elements_wt_percent'] == {'C': 50.0, 'H': 6.0, 'O': 44.0}
    assert daf['ash'] == {'value': 5.0, 'unit': '%'} and daf['total_kg'] == 1.0
    no_ash = legacy_feed_to_v2({'elements_wt_percent': {'C': 50.0, 'H': 6.0, 'O': 44.0}})
    assert 'ash' not in no_ash
    with pytest.raises(LegacyConversionError, match='unmodelled_wt_percent'):
        legacy_feed_to_v2({'elements_wt_percent': {'C': 50.0, 'H': 6.0, 'O': 34.0}, 'ash_wt_percent': 5.0,
                           'unmodelled_wt_percent': 5.0})
    with pytest.raises(LegacyConversionError, match='species data'):
        legacy_feed_to_v2({'species_moles': {'cellulose': 1.0}, 'water_g_per_g_dry': 4.0})


def test_legacy_species_feed():
    molar_mass = molar_mass_function(species_elements)
    v2 = legacy_feed_to_v2({'species_moles': {'cellulose': 2.0, 'methanol': 1.0}, 'water_g_per_g_dry': 4.0},
                           molar_mass)
    dry_g = 2 * molar_mass('cellulose') + molar_mass('methanol')
    assert v2['species'] == {'water': {'value': pytest.approx(4.0 * dry_g / 1000), 'unit': 'kg'},
                             'cellulose': {'value': 2.0, 'unit': 'mol'},
                             'methanol': {'value': 1.0, 'unit': 'mol'}}
    model = feed_v2_to_model(v2, molar_mass=molar_mass)
    assert model.species_moles == {'cellulose': 2.0, 'methanol': 1.0}
    assert model.water_g_per_g_dry == pytest.approx(4.0, rel=1e-14)
    deferred = feed_v2_to_model(v2)
    assert deferred.species_amounts == v2['species']
    with pytest.raises(ValueError, match='resolved'):
        deferred.total_water_g_per_g_dry()
    assert deferred.resolved(molar_mass) == model


def test_legacy_carbon_and_starts():
    data = legacy(constraints={'gas_carbon_min_fraction': 0.9, 'solid_carbon_fraction': 0.05,
                               'phase_carbon_fractions': {'organic_liquid': 0.01}},
                  solver={'method': 'ipopt', 'multistart': 4})
    v2 = legacy_to_v2(data)
    assert v2['format'] == 2 and list(v2)[0] == 'format'
    assert v2['constraints'] == {'carbon': {'gas': {'min': 0.9}, 'solid': {'exact': 0.05},
                                            'organic_liquid': {'exact': 0.01}}}
    # gas, aqueous and organic liquid have species (the cut-off at 390 °C is above 300 °C)
    assert v2['solver'] == {'method': 'ipopt', 'starts': 1 + 2 * 3 + 4}
    assert legacy_to_v2(legacy(solver={'multistart': 0}))['solver'] == {'starts': 1}
    hot = legacy(solver={'multistart': 2})
    hot['conditions']['temperature_C'] = [400, 450]  # aqueous removed at every temperature
    assert legacy_to_v2(hot)['solver'] == {'starts': 1 + 2 * 2 + 2}
    # the same numbers in a Case built from the old keys
    assert case_from_dict(data).starts == 11 and case_from_dict(hot).starts == 7
    assert legacy_to_v2(v2) == v2  # already format 2


def test_refine_is_merged_into_the_temperatures():
    data = legacy(sweep={'continuation': True, 'refine': {'start_C': 340, 'stop_C': 370, 'step_C': 10}})
    v2 = legacy_to_v2(data)
    assert v2['conditions']['temperature_C'] == [300, 340, 350, 360, 370, 400, 450]
    assert v2['sweep'] == {'continuation': True}
    assert case_from_dict(v2).temperatures() == case_from_dict(data).temperatures()
    only = legacy_to_v2(legacy(sweep={'refine': {'start_C': 340, 'stop_C': 370, 'step_C': 10}}))
    assert 'sweep' not in only


def test_legacy_carbon_bounds_mapping():
    assert legacy_carbon_bounds(gas_carbon_fraction=0.8) == ({'gas': {'exact': 0.8}}, 'char')
    assert legacy_carbon_bounds(gas_carbon_min_fraction=0.9, phase_carbon_fractions={'aqueous': 0.0}) == \
        ({'gas': {'min': 0.9}, 'aqueous': {'exact': 0.0}}, 'char')
    assert legacy_carbon_bounds(phase_carbon_fractions={'solid': 0.1}) == ({'solid': {'exact': 0.1}}, 'all_solids')
    with pytest.raises(ValueError, match='not both'):
        legacy_carbon_bounds(gas_carbon_fraction=0.8, gas_carbon_min_fraction=0.9)
    with pytest.raises(ValueError, match='Conflicting'):
        legacy_carbon_bounds(gas_carbon_fraction=0.8, phase_carbon_fractions={'gas': 0.7})


# ---------------------------------------------------------------- the Case
def test_case_moves_old_settings_into_the_new_ones():
    case = htg.load_case(TEMPLATE).replace(gas_carbon_fraction=0.8, phase_carbon_fractions={'aqueous': 0.0},
                                           multistart=2)
    assert case.carbon_bounds == {'gas': {'exact': 0.8}, 'aqueous': {'exact': 0.0}}
    assert case.gas_carbon_fraction is None and case.phase_carbon_fractions is None
    assert case.starts == 1 + 2 * 3 + 2 and case.multistart is None
    assert case.replace() == case
    with pytest.raises(ValueError, match='given twice'):
        case.replace(gas_carbon_min_fraction=0.5)


@pytest.mark.parametrize('bounds,message', [
    ({'gas': {'exact': 0.5, 'min': 0.2}}, 'Gas carbon: exact cannot be combined with min or max.'),
    ({'solid': {'min': 0.4, 'max': 0.2}}, 'Solid: the minimum is larger than the maximum.'),
    ({'aqueous': {'max': 1.5}}, 'Aqueous carbon: max must be a fraction between 0 and 1'),
    ({'organic_liquid': {'at_least': 0.1}}, 'Organic liquid carbon: use exact, or min and/or max'),
    ({'liquid': {'exact': 0.1}}, "Unknown phase 'liquid'"),
])
def test_carbon_bound_validation(bounds, message):
    with pytest.raises(ValueError, match=message.replace('(', r'\(')):
        check_carbon_bounds(bounds)
    with pytest.raises(ValueError, match=message.replace('(', r'\(')):
        htg.load_case(TEMPLATE).replace(carbon_bounds=bounds).validate()


def test_carbon_bounds_in_a_file(tmp_path):
    text = TEMPLATE.read_text().replace('# [constraints.carbon]\n# gas = { min = 0.9 }\n# solid = { max = 0.05 }',
                                        '[constraints.carbon]\ngas = { min = 0.9 }\nsolid = { max = 0.05 }')
    (tmp_path / 'c.toml').write_text(text)
    case = htg.load_case(tmp_path / 'c.toml')
    assert case.carbon_bounds == {'gas': {'min': 0.9}, 'solid': {'max': 0.05}}
    (tmp_path / 'bad.toml').write_text(text.replace('solid = { max = 0.05 }', 'solid = { min = 0.2, max = 0.05 }'))
    with pytest.raises(ValueError, match='Solid: the minimum is larger than the maximum.'):
        htg.load_case(tmp_path / 'bad.toml')


def test_starts_validation():
    case = htg.load_case(TEMPLATE)
    assert case.starts == 1 and case.distinct_composition_tol == 1e-4 and case.distinct_gibbs_tol == 1.0
    for bad in (0, 2.5, 'x'):
        with pytest.raises(ValueError, match='starts'):
            case.replace(starts=bad).validate()
    with pytest.raises(ValueError, match='either starts or the old multistart'):
        case.replace(starts=3, multistart=1)


def test_starts_from_multistart():
    assert starts_from_multistart(0, 3) == 1
    assert starts_from_multistart(4, 2) == 9
    phases = {'gas': ['water'], 'aqueous': ['water'], 'organic_liquid': [], 'solid': ['carbon']}
    assert fluid_phase_count(phases) == 2
    assert fluid_phase_count(phases, aqueous_enabled=False) == 1
    assert fluid_phase_count(phases, aqueous_cutoff_C=390, temperatures=[300, 400]) == 2
    assert fluid_phase_count(phases, aqueous_cutoff_C=390, temperatures=[390, 400]) == 1


def test_mixed_feed_keys_are_refused(tmp_path):
    text = TEMPLATE.read_text().replace('total_kg = 5.0', 'total_kg = 5.0\nwater_g_per_g_dry = 4.0')
    (tmp_path / 'c.toml').write_text(text)
    with pytest.raises(ValueError, match='mixes the new keys'):
        htg.load_case(tmp_path / 'c.toml')
    data = tomllib.loads(TEMPLATE.read_text())
    data['feed'] = {'elements_wt_percent': data['feed']['elements_wt_percent']}
    with pytest.raises(ValueError, match='total_kg and solids'):
        case_from_dict(data)
    del data['format']  # an old file: elements only, no water
    assert case_from_dict(data).feed.water_g_per_g_dry == 0


def test_species_feed_in_a_file(tmp_path):
    text = TEMPLATE.read_text()
    start, end = text.index('total_kg = 5.0'), text.index('# Or, instead')
    text = text[:start] + text[end:]
    text = text.replace('# [feed.species]\n# water = { value = 4.5, unit = "kg" }      # always present; "kg" or "mol"\n'
                        '# cellulose = { value = 0.5, unit = "kg" }',
                        '[feed.species]\nwater = { value = 4.5, unit = "kg" }\ncellulose = { value = 0.5, unit = "kg" }')
    text = text.replace('[feed]\n# The wet feed going into the reactor: solids plus water.\n', '')
    (tmp_path / 'c.toml').write_text(text)
    case = htg.load_case(tmp_path / 'c.toml')
    assert case.feed_v2 == {'species': {'water': {'value': 4.5, 'unit': 'kg'},
                                        'cellulose': {'value': 0.5, 'unit': 'kg'}}}
    resolved = case.feed.resolved(molar_mass_function(species_elements))
    assert resolved.water_g_per_g_dry == pytest.approx(9.0)
    assert resolved.species_moles['cellulose'] == pytest.approx(500 / molar_mass_function(species_elements)('cellulose'))
    case.save(tmp_path / 'again.toml')
    assert htg.load_case(tmp_path / 'again.toml').feed_v2 == case.feed_v2


# ---------------------------------------------------------------- feed arithmetic
FEED = {'total_kg': 10.0, 'solids': {'value': 15.0, 'unit': '%'}, 'ash': {'value': 8.45, 'unit': '%'},
        'elements_wt_percent': {'C': 30.4, 'H': 4.41, 'O': 53.58, 'N': 2.66, 'S': 0.51}}


def test_feed_summary_matches_the_case_tab_example():
    summary = feed_summary(FEED)
    assert summary['mode'] == 'elements'
    assert {e: round(v, 2) for e, v in summary['used'].items()} == \
        {'C': 33.20, 'H': 4.82, 'O': 58.52, 'N': 2.91, 'S': 0.56}
    assert round(summary['entered_sum'], 2) == 91.56
    assert round(summary['water_kg'], 3) == 8.500 and round(summary['water_mol'], 1) == 471.8
    assert round(summary['ash_kg'], 3) == 0.127 and round(summary['ash_free_kg'], 3) == 1.373
    assert round(summary['elements']['C']['kg'], 3) == 0.456
    assert round(summary['elements']['C']['mol'], 2) == 37.96
    assert [row['name'] for row in summary['rows']] == ['Total', 'Water', 'Solids', 'Ash', 'Ash-free solids',
                                                       'C', 'H', 'O', 'N', 'S']
    assert summary['rows'][0] == {'name': 'Total', 'kg': 10.0, 'mol': None}
    kg = feed_summary({**FEED, 'solids': {'value': 1.5, 'unit': 'kg'}, 'ash': {'value': 0.12675, 'unit': 'kg'}})
    assert kg['ash_free_kg'] == pytest.approx(summary['ash_free_kg'])


def test_feed_summary_species_mode():
    feed = {'species': {'water': {'value': 4.5, 'unit': 'kg'}, 'cellulose': {'value': 0.5, 'unit': 'kg'},
                        'methanol': {'value': 2.0, 'unit': 'mol'}}}
    summary = feed_summary(feed, species_elements=species_elements)
    M = molar_mass_function(species_elements)
    assert summary['mode'] == 'species'
    names = [row['name'] for row in summary['rows']]
    assert names[:5] == ['Total', 'Water', 'Solids (not water)', 'cellulose', 'methanol']
    assert set(names[5:]) == {'C', 'H', 'O'}
    dry = 0.5 + 2.0 * M('methanol') / 1000
    assert summary['solids_kg'] == pytest.approx(dry)
    assert summary['total_kg'] == pytest.approx(4.5 + dry)
    carbon = 6 * 500 / M('cellulose') + 2.0
    assert summary['elements']['C']['mol'] == pytest.approx(carbon)
    with pytest.raises(FeedError, match='species data'):
        feed_summary(feed)


@pytest.mark.parametrize('change,message', [
    ({'ash': {'value': 2.0, 'unit': 'kg'}}, "Ash can't be more than the solids."),
    ({'elements_wt_percent': {}}, 'Elements: give at least one element.'),
    ({'solids': {'value': 15.0, 'unit': 'g'}}, 'Solids: the unit must be'),
    ({'total_kg': 0}, 'Total mass must be positive.'),
    ({'solids': {'value': 12.0, 'unit': 'kg'}}, "Solids can't be more than the total mass."),
])
def test_feed_errors(change, message):
    with pytest.raises(FeedError, match=message):
        feed_summary({**FEED, **change})
    with pytest.raises(FeedError, match=message):
        feed_v2_to_model({**FEED, **change})


def test_feed_v2_to_model():
    model = feed_v2_to_model(FEED)
    assert model.basis == 'dry_ash_free' and model.ash_wt_percent == 8.45 and model.moisture_wt_percent == 0
    assert model.water_g_per_g_dry == pytest.approx(8.5 / 1.5)
    assert sum(model.elements_wt_percent.values()) == pytest.approx(100)
    assert model.elements_wt_percent['C'] == pytest.approx(30.4 / 91.56 * 100)
    dry = model.dry_element_wt_percent(None)
    assert sum(dry.values()) == pytest.approx(100 - 8.45)
    in_kg = feed_v2_to_model({**FEED, 'ash': {'value': 0.12675, 'unit': 'kg'}})
    assert in_kg.ash_wt_percent == pytest.approx(8.45)
    no_ash = feed_v2_to_model({k: v for k, v in FEED.items() if k != 'ash'})
    assert no_ash.ash_wt_percent is None
    s, a = feed_fractions(model)
    assert s == pytest.approx(0.15) and a == pytest.approx(0.0845)
    salts = feed_v2_to_model({**FEED, 'salts': {'KCl': {'mol_per_g_dry_feed': 1e-5, 'elements': {'K': 1, 'Cl': 1}}}})
    assert salts.salts['KCl'].mol_per_g_dry_feed == 1e-5


def test_feed_fractions_of_old_feeds():
    old = htg.load_case(LEGACY / 'yakaboylu2013_fig7_8.toml').feed
    s, a = feed_fractions(old)
    assert s == pytest.approx(1 / 4.999935304) and a == pytest.approx(0.0844688764)


# ---------------------------------------------------------------- reproducing saved results
SAVED = Path(os.environ.get('HTG_SAVED_RESULTS', ROOT / 'results'))


@pytest.mark.slow
@pytest.mark.parametrize('name', [p.stem for p in sorted(CASES.glob('*.toml'))])
def test_converted_cases_reproduce_saved_results(name):
    """Run each converted example case and compare with results/<case>/summary.csv.

    Set HTG_SAVED_RESULTS to a results folder made before format 2 (skipped if missing).
    """
    folder = SAVED / name
    if not (folder / 'summary.csv').exists():
        pytest.skip(f'no saved results in {folder}')
    from htg.runner import WorkerPool, default_workers
    case = htg.load_case(CASES / f'{name}.toml')
    saved = pd.read_csv(folder / 'summary.csv')
    with WorkerPool(min(4, default_workers(), len(case.chains())), {'data_dir': str(ROOT / 'data')},
                    warm_up=False) as pool:
        run = pool.submit(case)
        assert run.wait(3 * 3600)
    summary = run.sweep_result().summary().set_index(['T_C', 'P_bar'])
    for _, row in saved.iterrows():
        new = summary.loc[(row.T_C, row.P_bar)]
        assert new.status == row.status, (row.T_C, new.get('error'))
        if row.status == 'solved':
            assert new.gibbs_J_per_mol_C == pytest.approx(row.gibbs_J_per_mol_C, rel=1e-9), row.T_C
