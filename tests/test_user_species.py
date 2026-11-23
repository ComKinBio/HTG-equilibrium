"""WP8: user species database, Species tab catalogue and filters, import from file.

None of these tests read the real user folder: every model gets ``user_data``
pointing into tmp_path (or None), and HTG_USER_DATA is redirected.
"""
import tomllib
from pathlib import Path

import pandas as pd
import pytest

import htg
from htg.cli import main
from htg.formulas import formula_summary, parse_formula
from htg.runner import WorkerPool
from htg.species_catalogue import (DATA_SOURCES, PLAIN_REASONS, build_catalogue, count_text,
                                   filter_catalogue)
from htg.species_import import NEW, REPLACES, SKIPPED, apply_import, detect_format, preview_import
from htg.user_species import (EPCSAFT_COLUMNS, PCSAFT_COLUMNS, RUN_RECORD_FILE, SpeciesSpec,
                              UserDatabase, check_species, format_nasa_record, nasa_match_message,
                              parse_nasa_record, spec_from_model, user_species_table)

DATA = Path(htg.__file__).resolve().parents[2] / 'data'
TEMPLATE = Path(htg.__file__).parent / 'templates' / 'case_template.toml'
T = 623.15


@pytest.fixture(autouse=True)
def no_real_user_folder(tmp_path, monkeypatch):
    monkeypatch.setenv('HTG_USER_DATA', str(tmp_path / 'env-user-folder'))


@pytest.fixture(scope='module')
def shipped():
    return htg.Model.load(user_data=None)


@pytest.fixture(scope='module')
def catalogue(shipped):
    return build_catalogue(shipped)


def methane_m(model):
    return model.registry.resolve_parameters(model.databases, T, ['methane'])['methane']['m']


def nasa_coefficients(model, record):
    r = model.databases.nasa[record]
    return {'T_low': r['T_low'], 'T_mid': r['T_mid'], 'T_high': r['T_high'],
            'high': list(r['high']), 'low': list(r['low'])}


# ---------------------------------------------------------------- formulas and NASA records
def test_formula_summary():
    assert formula_summary('CH3OH') == 'C 1 · H 4 · O 1 · 32.04 g/mol'
    assert parse_formula('CO3^2-') == ({'C': 1, 'O': 3}, -2)
    assert parse_formula('Mg2+', charge=2) == ({'Mg': 1}, 2)
    assert parse_formula('CaSO4*2H2O')[0] == {'Ca': 1, 'S': 1, 'O': 6, 'H': 4}
    with pytest.raises(ValueError, match='Unknown element'):
        parse_formula('Xq2')


def test_nasa_record_round_trip(shipped):
    c = nasa_coefficients(shipped, 'CH4')
    text = format_nasa_record('MYCH4', {'C': 1, 'H': 4}, 'G', c['T_low'], c['T_mid'], c['T_high'],
                              c['high'], c['low'])
    assert all(len(line) == 80 for line in text.splitlines())
    parsed = parse_nasa_record(text)
    assert parsed['name'] == 'MYCH4' and parsed['elements'] == {'C': 1.0, 'H': 4.0}
    assert parsed['high'] == pytest.approx(c['high']) and parsed['low'] == pytest.approx(c['low'])
    # Pasted without the column-80 markers and trailing spaces.
    lines = text.splitlines()
    loose = '\n'.join([lines[0][:73]] + [line[:75].rstrip() for line in lines[1:]])
    assert parse_nasa_record(loose)['low'] == pytest.approx(c['low'])
    with pytest.raises(ValueError, match='4 lines'):
        parse_nasa_record('\n'.join(loose.splitlines()[:3]))


def test_shipped_columns_match():
    assert list(pd.read_csv(DATA / 'PROJECT_PCSAFT.csv').columns) == PCSAFT_COLUMNS
    assert list(pd.read_csv(DATA / 'PROJECT_EPCSAFT.csv').columns) == EPCSAFT_COLUMNS


# ---------------------------------------------------------------- catalogue and filters
def test_catalogue_columns_and_status(catalogue):
    assert len(catalogue) == catalogue.Species.nunique()
    methane = catalogue.set_index('Species').loc['methane']
    assert methane.Usable and methane.Status == 'Usable'
    assert methane['Can be in'] == ['gas', 'aqueous', 'organic_liquid']
    assert {'NASA-7', 'PC-SAFT'} <= set(methane.Data) and not methane.Yours
    carbon = catalogue.set_index('Species').loc['carbon']
    assert carbon['Can be in'] == ['solid'] and carbon.Usable
    assert catalogue.set_index('Species').loc['Na+']['Can be in'] == ['aqueous']
    # Only plain-language reasons, never the raw registry texts.
    raw = set(PLAIN_REASONS) - {'NASA solid'}
    assert not set(catalogue.Status) & raw
    assert catalogue.Status.str.startswith('No NASA-7 record with formula ').any()
    assert (catalogue.Status == 'No ePC-SAFT ion parameters').any()


def test_filter_contains_elements(catalogue):
    shown = filter_catalogue(catalogue, elements='C, H', show='All')
    names = set(shown.Species)
    assert {'methane', 'methanol', 'benzene'} <= names
    assert 'carbon dioxide' not in names and 'water' not in names
    assert all({'C', 'H'} <= set(e) for e in shown.Elements)


def test_filter_exact_match(catalogue):
    shown = filter_catalogue(catalogue, elements='c,h', exact=True, show='All')
    names = set(shown.Species)
    assert {'methane', 'benzene'} <= names and 'methanol' not in names
    assert all(set(e) == {'C', 'H'} for e in shown.Elements)
    cho = set(filter_catalogue(catalogue, elements='C H O', exact=True).Species)
    assert 'methanol' in cho and 'methane' not in cho
    with pytest.raises(ValueError, match='Unknown element'):
        filter_catalogue(catalogue, elements='C, Qq')


def test_filter_can_be_in(catalogue):
    solids = filter_catalogue(catalogue, can_be_in=['solid'])
    assert len(solids) and all('solid' in t for t in solids['Can be in'])
    assert 'carbon' in set(solids.Species) and 'methane' not in set(solids.Species)
    both = filter_catalogue(catalogue, can_be_in=['Gas', 'Solid'])
    assert {'carbon', 'methane'} <= set(both.Species)
    aqueous = set(filter_catalogue(catalogue, can_be_in=['aqueous']).Species)
    assert {'Na+', 'methane'} <= aqueous and 'carbon' not in aqueous


def test_filter_show_source_text_and_count(catalogue):
    usable = filter_catalogue(catalogue)
    unusable = filter_catalogue(catalogue, show='Not usable only')
    every = filter_catalogue(catalogue, show='All')
    assert usable.Usable.all() and not unusable.Usable.any()
    assert len(usable) + len(unusable) == len(every) == len(catalogue)
    supcrt = filter_catalogue(catalogue, source='SUPCRT92', show='All')
    assert all('SUPCRT92' in d for d in supcrt.Data) and set(DATA_SOURCES) >= {'SUPCRT92'}
    meth = filter_catalogue(catalogue, text='METHAN')
    assert 'methane' in set(meth.Species) and 'methanol' in set(meth.Species)
    assert count_text(meth, catalogue) == f'{len(meth)} of {len(catalogue)} species'
    assert count_text(6, 412) == '6 of 412 species'


# ---------------------------------------------------------------- user database
def test_new_molecule_with_complete_data_is_usable(shipped, tmp_path):
    folder = tmp_path / 'user'
    spec = SpeciesSpec(kind='molecule', name='Testanol', formula='C9H12O3', thermo='nasa',
                       nasa=nasa_coefficients(shipped, 'CH4'), model='PC-SAFT', m=3.1, sigma=3.6,
                       epsilon=250.0, association_scheme='2B', epsilonAB=2000.0, kappaAB=0.02,
                       kij=[('methane', 0.01)], authors='A. Author', year='2024', doi='10.1/x')
    check = check_species(spec, shipped)
    assert check.usable, check.items
    assert [i.text for i in check.items] == ['Formula understood', 'Thermochemistry complete',
                                             'Equation of state parameters complete',
                                             'Name not used by another species']
    assert check.can_be_in_text == 'Can be in: gas, aqueous, organic liquid'
    name = UserDatabase.load(folder).save_species(spec, shipped)
    assert name == 'testanol'
    model = htg.Model.load(user_data=folder)
    row = build_catalogue(model).set_index('Species').loc['testanol']
    assert row.Usable and row.Yours and not row.Overrides and 'PC-SAFT' in row.Data
    params = model.registry.resolve_parameters(model.databases, T, ['testanol'])['testanol']
    assert params['m'] == 3.1 and params['Source'] == 'user_pcsaft.csv'
    assert (('C9H12O3', 'CH4'), 0.01) in model.registry.binary_interactions
    assert 'testanol' not in set(build_catalogue(shipped).Species)
    table = user_species_table(model)
    assert table.set_index('Species').loc['testanol'].Origin == 'yours'
    # Editing it again starts from the stored values.
    again = spec_from_model(model, 'testanol')
    assert again.m == 3.1 and again.authors == 'A. Author' and again.thermo == 'nasa_record'


def test_override_methane_m_and_switch_off(shipped, tmp_path):
    folder = tmp_path / 'user'
    original = methane_m(shipped)
    spec = spec_from_model(shipped, 'methane')
    spec.m = original * 1.5
    assert check_species(spec, shipped, editing='methane').usable
    UserDatabase.load(folder).save_species(spec, shipped)
    model = htg.Model.load(user_data=folder)
    assert methane_m(model) == pytest.approx(original * 1.5)
    assert model.user_info['methane']['overrides']
    assert build_catalogue(model).set_index('Species').loc['methane'].Overrides
    assert model.load_args['user_data'] == str(folder)  # worker processes load the same layer
    # Without the user layer (user_data=None, i.e. --no-user-species) the shipped value is back.
    assert methane_m(htg.Model.load(user_data=None)) == original
    # Reset to original.
    UserDatabase.load(folder).reset_species('methane')
    assert methane_m(htg.Model.load(user_data=folder)) == original


def test_cli_user_species_flags(shipped, tmp_path, capsys):
    folder = tmp_path / 'user'
    spec = spec_from_model(shipped, 'methane')
    spec.name, spec.formula = 'clitestane', 'C11H24O'
    spec.thermo, spec.nasa = 'nasa', nasa_coefficients(shipped, 'CH4')
    spec.nasa_record = ''
    UserDatabase.load(folder).save_species(spec, shipped)
    assert main(['species', '--user-data', str(folder), '--search', 'clitestane']) == 0
    assert 'clitestane' in capsys.readouterr().out
    assert main(['species', '--user-data', str(folder), '--no-user-species',
                 '--search', 'clitestane']) == 0
    assert '0 species' in capsys.readouterr().out


def test_ion_and_solids(shipped, tmp_path):
    folder = tmp_path / 'user'
    db = UserDatabase.load(folder)
    ion = SpeciesSpec(kind='ion', name='Ag+', formula='Ag+', charge=1, thermo='hkf_record',
                      supcrt_record='Ag+', m=1.0, sigma=2.9, epsilon=200.0)
    assert check_species(ion, shipped, editing='Ag+').usable  # completes the SUPCRT-only Ag+
    db.save_species(ion, shipped)
    mineral = SpeciesSpec(kind='solid', name='Testite', formula='CaCO3', thermo='maier_kelley',
                          maier_kelley={'Gf': -269880.0, 'Hf': -288552.0, 'Sr': 22.15, 'Vr': 36.934,
                                        'a': 24.98, 'b': 0.00524, 'c': -620000.0, 'Tmax': 1200.0})
    assert check_species(mineral, shipped).usable
    db.save_species(mineral, shipped)
    nasa_solid = SpeciesSpec(kind='solid', name='mychar', formula='C', thermo='nasa',
                             nasa=nasa_coefficients(shipped, 'CHAR'))
    assert check_species(nasa_solid, shipped).usable
    db.save_species(nasa_solid, shipped)
    model = htg.Model.load(user_data=folder)
    table = build_catalogue(model).set_index('Species')
    assert table.loc['Ag+'].Usable and table.loc['Ag+']['Can be in'] == ['aqueous']
    assert table.loc['Testite'].Usable and table.loc['Testite']['Can be in'] == ['solid']
    assert table.loc['mychar'].Usable and table.loc['mychar'].Yours
    assert model.registry.ion['Ag+']['charge'] == 1
    assert model.registry.solid['Testite']['supcrt_mineral']['Tmax'] == 1200.0
    # Remove from my database.
    UserDatabase.load(folder).remove_species('Testite')
    assert 'Testite' not in htg.Model.load(user_data=folder).registry.solid


def test_incomplete_species_is_saved_but_not_usable(shipped, tmp_path):
    folder = tmp_path / 'user'
    spec = SpeciesSpec(kind='molecule', name='Methane', formula='C2H6O', thermo='nasa',
                       nasa={'T_low': 300, 'T_mid': 1000}, m=None, sigma=3.0, epsilon=None)
    check = check_species(spec, shipped)
    texts = [i.text for i in check.items if not i.ok]
    assert not check.usable
    assert any(t.startswith('Missing NASA-7 values: T high') for t in texts)
    assert 'Missing equation of state parameters: m, ε/k' in texts
    assert any(t.startswith('Name already used by methane') for t in texts)
    spec.name = 'halfdone'
    UserDatabase.load(folder).save_species(spec, shipped)
    row = build_catalogue(htg.Model.load(user_data=folder)).set_index('Species').loc['halfdone']
    assert not row.Usable and row.Yours
    reopened = spec_from_model(htg.Model.load(user_data=folder), 'halfdone')
    assert reopened.thermo == 'nasa' and reopened.nasa['T_low'] == 300


def test_hidden_shipped_species(tmp_path):
    folder = tmp_path / 'user'
    UserDatabase.load(folder).hide_species('ethane')
    model = htg.Model.load(user_data=folder)
    assert 'ethane' not in model.registry.molecular
    assert 'ethane' not in set(build_catalogue(model).Species)
    UserDatabase.load(folder).hide_species('ethane', hidden=False)
    assert 'ethane' in htg.Model.load(user_data=folder).registry.molecular


def test_nasa_match_message(shipped):
    assert nasa_match_message(shipped, 'CH4') == 'Only record with formula CH4'
    assert nasa_match_message(shipped, 'C9H12O3').startswith('No record with formula')


# ---------------------------------------------------------------- import
def _pcsaft_csv():
    ship = pd.read_csv(DATA / 'PROJECT_PCSAFT.csv')
    base = ship[(ship.Species == 'methane') & (ship.Model == 'PC-SAFT')].head(1)
    new = base.assign(Species='importane', Formula='C13H28O', m=4.2)
    existing = base.assign(m=1.234)
    broken = base.assign(Species='brokane', m='not a number')
    return pd.concat([new, existing, broken]).to_csv(index=False)


def test_import_pcsaft_table(shipped, tmp_path):
    folder = tmp_path / 'user'
    model = htg.Model.load(user_data=folder)
    text = _pcsaft_csv()
    assert detect_format(text) == 'pcsaft'
    rows = preview_import(text, model)
    assert [r.action for r in rows] == [NEW, REPLACES, SKIPPED]
    new, existing, broken = rows
    assert new.selected and not existing.selected and not broken.selectable
    assert 'm' in existing.differs and 'differs: ' in existing.after
    assert broken.reason == 'm is not a positive number'
    assert new.after.startswith('No NASA-7 record with formula C13H28O')
    assert apply_import([r for r in rows if r.selected], model.user_data) == 1
    after = htg.Model.load(user_data=folder)
    table = build_catalogue(after).set_index('Species')
    assert table.loc['importane'].Yours and 'brokane' not in table.index
    assert methane_m(after) == methane_m(shipped)  # the replacement was not ticked


def test_detect_other_formats(shipped):
    lines = (DATA / 'Thermdat_WasteMech.txt').read_text().splitlines(keepends=True)
    assert detect_format(''.join(lines[1:9])) == 'nasa'
    rows = preview_import(''.join(lines[1:5]), shipped)
    assert rows[0].action == SKIPPED and rows[0].reason == 'same values as the existing record'
    xml = ('<Database><Species><Name>Testion+</Name><Formula>Ag+</Formula><Elements>Ag(1)</Elements>'
           '<Type>Aqueous</Type><Charge>1</Charge><Thermo><HKF><Gf>1</Gf></HKF></Thermo></Species>'
           '</Database>')
    assert detect_format(xml) == 'supcrt'
    assert preview_import(xml, shipped)[0].reason == 'incomplete HKF data'


# ---------------------------------------------------------------- saved runs
def fake_worker(worker_id, inbox, outbox, model_args, warm_up):
    outbox.put(('state', worker_id, 'idle', ''))
    while True:
        message = inbox.get()
        if message is None:
            break
        task_id, case, chain = message
        for T_, P in chain:
            outbox.put(('point', worker_id, task_id, T_, P, 'fake', 'test'))
        outbox.put(('task_done', worker_id, task_id))


def test_saved_runs_record_user_species(shipped, tmp_path):
    folder = tmp_path / 'user'
    spec = spec_from_model(shipped, 'methane')
    spec.m = 1.5
    UserDatabase.load(folder).save_species(spec, shipped)
    case = htg.load_case(TEMPLATE).replace(name='t', temperature_C=[300], pressure_bar=240)
    assert 'methane' in {s for names in case.phases.values() for s in names}
    with WorkerPool(1, {'user_data': str(folder)}, worker_target=fake_worker) as pool:
        run = pool.submit(case, output_dir=tmp_path / 'run')
        assert run.wait(20)
    out = tmp_path / 'run'
    with open(out / RUN_RECORD_FILE, 'rb') as f:
        record = tomllib.load(f)
    assert record['database'] == str(folder.resolve())
    methane = record['species']['methane']
    assert methane['epcsaft'][0]['Value'] == '1.5' or float(methane['epcsaft'][0]['Value']) == 1.5
    assert f'# User species database: {folder.resolve()}' in (out / 'case.toml').read_text()
    htg.load_case(out / 'case.toml')  # the note keeps the case file loadable
    # Without the user layer, no record and no note.
    with WorkerPool(1, {'user_data': None}, worker_target=fake_worker) as pool:
        plain = pool.submit(case, output_dir=tmp_path / 'plain')
        assert plain.wait(20)
    assert not (tmp_path / 'plain' / RUN_RECORD_FILE).exists()
    assert 'User species' not in (tmp_path / 'plain' / 'case.toml').read_text()


@pytest.mark.slow
def test_sweep_with_user_species_solves_and_records(shipped, tmp_path):
    """A user copy of methane (same values) gives the shipped result; the saved run records it."""
    folder = tmp_path / 'user'
    UserDatabase.load(folder).save_species(spec_from_model(shipped, 'methane'), shipped)
    model = htg.Model.load(user_data=folder)
    case = htg.load_case(TEMPLATE).replace(name='t', temperature_C=[400], pressure_bar=240)
    sweep = model.sweep(case, progress=None)
    reference = shipped.sweep(case, progress=None)
    assert sweep.results and reference.results
    assert sweep.results[0].gibbs_J_per_mol_C == pytest.approx(
        reference.results[0].gibbs_J_per_mol_C, rel=1e-8)
    sweep.save(tmp_path / 'run')
    with open(tmp_path / 'run' / RUN_RECORD_FILE, 'rb') as f:
        assert 'methane' in tomllib.load(f)['species']
