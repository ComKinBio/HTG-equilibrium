"""Carbon bounds, the aqueous cut-off rule, starting points and distinct results.

Most tests need no Julia: they use the pure helper functions, a stand-in
registry, or a small stand-in Gibbs problem with a scripted optimizer backend.
"""
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import htg
from htg import optimizer as opt
from htg.case import carbon_bound_messages, effective_carbon_bounds
from htg.equilibrium import apply_phase_rules, carbon_bound_rows
from htg.results import (Alternative, Result, SweepResult, apply_choices, read_choices,
                         write_choices)
from htg.gui.runs import live_tables, load_saved

TEMPLATE = Path(htg.__file__).parent / 'templates' / 'case_template.toml'
PAIRS = [('gas', 'water'), ('gas', 'methane'), ('gas', 'carbon dioxide'), ('aqueous', 'water'),
         ('aqueous', 'methane'), ('organic_liquid', 'benzene'), ('solid', 'carbon'), ('solid', 'Magnesite')]
CARBON = np.array([0, 1, 1, 0, 1, 6, 1, 1], dtype=float)


def case(**changes):
    return htg.load_case(TEMPLATE).replace(**changes)


def mask(*phases):
    return CARBON * np.array([p in phases for p, _ in PAIRS], dtype=float)


# ---------------------------------------------------------------- carbon bound rows
def test_rows_for_each_kind_of_bound():
    exact, bounds = carbon_bound_rows({'gas': {'min': 0.6, 'max': 0.9}, 'solid': {'max': 0.05},
                                       'aqueous': {'exact': 0.03}, 'organic_liquid': {'min': 0.01}},
                                      PAIRS, CARBON, 'char')
    assert [(p, f) for p, _, f in exact] == [('aqueous', 0.03)]
    np.testing.assert_array_equal(exact[0][1], mask('aqueous'))
    assert [(p, k, f) for p, k, _, f in bounds] == [('gas', 'min', 0.6), ('gas', 'max', 0.9),
                                                     ('solid', 'max', 0.05), ('organic_liquid', 'min', 0.01)]
    np.testing.assert_array_equal(bounds[0][2], mask('gas'))
    np.testing.assert_array_equal(bounds[3][2], mask('organic_liquid'))
    char = np.zeros(len(PAIRS))
    char[PAIRS.index(('solid', 'carbon'))] = 1.0
    np.testing.assert_array_equal(bounds[2][2], char)  # char only


def test_solid_rows_follow_the_target():
    _, bounds = carbon_bound_rows({'solid': {'min': 0.1}}, PAIRS, CARBON, 'all_solids')
    np.testing.assert_array_equal(bounds[0][2], mask('solid'))
    with pytest.raises(ValueError, match="Select solid 'carbon'"):
        carbon_bound_rows({'solid': {'exact': 0.1}}, PAIRS[:-2], CARBON[:-2], 'char')


def test_legacy_gas_minimum_gives_the_old_row():
    """The old code used one inequality: gas carbon >= fraction * feed carbon."""
    old = case(gas_carbon_min_fraction=0.9)
    assert old.carbon_bounds == {'gas': {'min': 0.9}}
    exact, bounds = carbon_bound_rows(old.carbon_bounds, PAIRS, CARBON, old.solid_carbon_target)
    assert exact == [] and len(bounds) == 1
    phase, kind, row, fraction = bounds[0]
    gas_carbon = CARBON * np.array([p == 'gas' for p, _ in PAIRS], dtype=float)  # the old gas_min_row
    assert (phase, kind, fraction) == ('gas', 'min', 0.9)
    np.testing.assert_array_equal(row, gas_carbon)
    exact, _ = carbon_bound_rows(case(gas_carbon_fraction=0.8, solid_carbon_fraction=0.1,
                                      phase_carbon_fractions={'organic_liquid': 0.05}).carbon_bounds,
                                 PAIRS, CARBON, 'char')
    assert [(p, f) for p, _, f in exact] == [('gas', 0.8), ('solid', 0.1), ('organic_liquid', 0.05)]


# ---------------------------------------------------------------- preflight messages
@pytest.mark.parametrize('bounds,message', [
    ({'gas': {'min': 0.7}, 'solid': {'exact': 0.4}}, 'The carbon shares add up to more than 100 %.'),
    ({'gas': {'exact': 0.5}, 'aqueous': {'exact': 0.2}, 'organic_liquid': {'exact': 0.1},
      'solid': {'exact': 0.1}}, 'With every phase fixed, the carbon shares must add up to 100 %.'),
    ({'gas': {'min': 0.6, 'max': 0.5}}, 'Gas: the minimum is larger than the maximum.'),
])
def test_preflight_messages(bounds, message):
    c = case(carbon_bounds=bounds, solid_carbon_target='all_solids')
    assert carbon_bound_messages(c) == [message]


def test_preflight_bound_on_empty_phase():
    c = case(carbon_bounds={'organic_liquid': {'max': 0.1}},
             phases={**case().phases, 'organic_liquid': []})
    assert carbon_bound_messages(c) == ['Organic liquid has a carbon bound but no species.']
    fine = case(carbon_bounds={'gas': {'exact': 0.5}, 'aqueous': {'exact': 0.2},
                               'organic_liquid': {'exact': 0.1}, 'solid': {'exact': 0.2}})
    assert carbon_bound_messages(fine) == []
    # with "char only", other solids may hold carbon, so the shares need not add up
    partial = fine.replace(carbon_bounds={**fine.carbon_bounds, 'solid': {'exact': 0.1}},
                           phases={**fine.phases, 'solid': ['carbon', 'Magnesite']})
    assert carbon_bound_messages(partial) == []
    assert carbon_bound_messages(partial, carbon_solids={'Magnesite': False}) == [
        'With every phase fixed, the carbon shares must add up to 100 %.']


# ---------------------------------------------------------------- aqueous cut-off
def test_cutoff_drops_the_aqueous_bound():
    registry = SimpleNamespace(ion={})
    c = case(carbon_bounds={'aqueous': {'exact': 0.03}, 'gas': {'min': 0.5}}, aqueous_cutoff_C=390)
    for T, dropped in ((350.0, []), (390.0, ['aqueous carbon']), (450.0, ['aqueous carbon'])):
        report = apply_phase_rules(c.at(T), registry)[-1]
        assert report['dropped_constraints'] == dropped
    assert effective_carbon_bounds(c.carbon_bounds, True) == ({'gas': {'min': 0.5}}, ['aqueous carbon'])
    assert effective_carbon_bounds({'aqueous': {'max': 0.1}}, True) == ({}, [])  # nothing to note
    assert effective_carbon_bounds(c.carbon_bounds, False) == (c.carbon_bounds, [])


# ---------------------------------------------------------------- grouping
def test_grouping_by_both_tolerances():
    x = np.array([1.0, 0.0, 0.5])
    runs = [('initial', x, -100.0), ('max gas', x + 5e-5, -100.4), ('random 0', x + [0, 0, 1e-3], -100.2),
            ('random 1', x, -98.0), ('min gas', x + 2e-5, -100.39)]
    groups = opt.group_solutions(runs, composition_tol=1e-4, gibbs_tol=1.0)
    # lowest G first: 'max gas' (-100.4); 'min gas', 'initial' are within both limits of it;
    # 'random 0' differs in composition; 'random 1' differs in G only
    assert [g['label'] for g in groups] == ['A', 'B', 'C']
    assert [[runs[i][0] for i in g['members']] for g in groups] == [
        ['max gas', 'min gas', 'initial'], ['random 0'], ['random 1']]
    assert [g['G'] for g in groups] == [-100.4, -100.2, -98.0]
    loose = opt.group_solutions(runs, composition_tol=1e-2, gibbs_tol=5.0)
    assert len(loose) == 1 and len(loose[0]['members']) == 5
    tight_g = opt.group_solutions(runs, composition_tol=1e-2, gibbs_tol=0.05)
    assert [[runs[i][0] for i in g['members']] for g in tight_g] == [
        ['max gas', 'min gas'], ['random 0'], ['initial'], ['random 1']]


def test_identical_runs_collapse():
    x = np.array([0.2, 0.3])
    groups = opt.group_solutions([('a', x, -5.0), ('b', x.copy(), -5.0), ('c', x.copy(), -5.0)])
    assert len(groups) == 1 and groups[0]['members'] == [0, 1, 2] and groups[0]['label'] == 'A'
    assert opt.result_label(0) == 'A' and opt.result_label(25) == 'Z' and opt.result_label(26) == 'AA'


# ---------------------------------------------------------------- the optimizer with a stand-in problem
class LinearPhase:
    """G = c . n (enough for the bookkeeping of starts and results)."""
    species = ['a', 'b']

    def __init__(self, c):
        self.c = np.asarray(c, dtype=float)

    def gibbs_and_mu(self, n):
        return float(self.c @ n), self.c.copy()

    def gibbs(self, n):
        return float(self.c @ n)

    def hessian(self, n):
        return np.zeros((2, 2))


def stand_in_problem():
    phase = LinearPhase([0.0, 10.0])
    return opt.GibbsProblem(
        phases={'gas': phase}, pairs=[('gas', 'a'), ('gas', 'b')], free=[0, 1],
        phase_indices={'gas': np.array([0, 1])}, ions=set(), T=600.0, A_eq=np.array([[1.0, 1.0]]),
        b_eq=np.array([1.0]), ineq_rows=None, ineq_targets=None, lower=[0.0, 0.0],
        x_start=np.array([0.5, 0.5]), G0=5.0, scale=1.0, model_total_G=lambda n: phase.gibbs(n))


def scripted(answers):
    """A backend returning the scripted solutions in turn (None: the start fails)."""
    calls = []

    def backend(problem, x0, options):
        x = answers[len(calls)]
        calls.append(x0)
        if x is None:
            return opt.OptimizeResult(x=x0, success=False, message='failed', nit=1)
        return opt.OptimizeResult(x=np.array(x, dtype=float), success=True, message='ok', nit=3)
    return backend, calls


def test_every_start_is_kept_and_grouped(monkeypatch):
    backend, calls = scripted([[0.99995, 0.00005], [1.0, 0.0], None, [0.0, 1.0], [0.3, 0.7]])
    monkeypatch.setitem(opt.BACKENDS, 'scripted', backend)
    problem = stand_in_problem()
    result, diagnostics = opt.optimize(problem, {'method': 'scripted', 'starts': 5, 'seed': 1})
    assert len(calls) == 5
    labels = [run['start'] for run in diagnostics['starts']]
    assert labels[:3] == ['initial', 'max gas', 'min gas'] and labels[3:] == ['random 0', 'random 1']
    groups = diagnostics['groups']
    assert [g['label'] for g in groups] == ['A', 'B', 'C']
    assert groups[0]['start_labels'] == ['max gas', 'initial'] and groups[0]['n_starts'] == 2
    assert groups[1]['G'] == pytest.approx(7.0) and groups[2]['G'] == pytest.approx(10.0)
    np.testing.assert_allclose(result.x, [1.0, 0.0])  # the main result is A
    assert 'tangent_plane_RT' in diagnostics  # stability test on by default with several starts


def test_one_start_is_the_old_behaviour(monkeypatch):
    backend, calls = scripted([[1.0, 0.0]])
    monkeypatch.setitem(opt.BACKENDS, 'scripted', backend)
    result, diagnostics = opt.optimize(stand_in_problem(), {'method': 'scripted'})
    assert len(calls) == 1 and [run['start'] for run in diagnostics['starts']] == ['initial']
    assert [g['label'] for g in diagnostics['groups']] == ['A']
    assert 'tangent_plane_RT' not in diagnostics


def test_candidate_starts_pad_with_random_mixes():
    problem = stand_in_problem()
    assert [label for label, _ in opt.candidate_starts(problem, 1)] == ['initial']
    assert [label for label, _ in opt.candidate_starts(problem, 2)] == ['initial', 'max gas']
    starts = opt.candidate_starts(problem, 6, seed=3)
    assert [label for label, _ in starts] == ['initial', 'max gas', 'min gas', 'random 0', 'random 1',
                                              'random 2']
    old = opt.starting_points(problem, n_random=3, seed=3)  # the old multistart = 3
    assert [label for label, _ in old] == [label for label, _ in starts]
    for (_, a), (_, b) in zip(old, starts):
        np.testing.assert_array_equal(a, b)
    for _, x in starts:
        assert x.sum() == pytest.approx(1.0) and (x >= -1e-12).all()


# ---------------------------------------------------------------- saved results and choices
def fake_result(T, alternatives=None, dropped=()):
    """A Result with two species, without solving."""
    amounts = {('gas', 'methane'): 0.6, ('gas', 'carbon dioxide'): 0.4, ('solid', 'carbon'): 0.0}
    species_data = {('gas', 'methane'): {'formula': 'CH4', 'elements': {'C': 1, 'H': 4}, 'charge': 0},
                    ('gas', 'carbon dioxide'): {'formula': 'CO2', 'elements': {'C': 1, 'O': 2}, 'charge': 0},
                    ('solid', 'carbon'): {'formula': 'C', 'elements': {'C': 1}, 'charge': 0}}
    result = Result(
        temperature_C=T, pressure_bar=240.0, amounts=amounts, initial_amounts=amounts, gibbs_J_per_mol_C=-100.0,
        phase_gibbs_J_per_mol_C={'gas': -100.0, 'aqueous': 0.0, 'organic_liquid': 0.0, 'solid': 0.0},
        balance_residual=0.0, charge_residual=0.0, iterations=3, eos_trial_errors=[],
        feed_C_mol_per_g_dry=0.03, elements=['C', 'H', 'O'], feed_atoms=np.array([1.0, 2.4, 0.8]),
        product_atoms=np.array([1.0, 2.4, 0.8]), gas_carbon_fraction=1.0, carbon_bound_violation=0.0,
        phase_rules={'aqueous_disabled': False, 'dropped_constraints': list(dropped)},
        optimizer={'method': 'ipopt', 'descent_RT': 0.0}, atomic_weights={'C': 12.011, 'H': 1.008, 'O': 15.999},
        species_data=species_data)
    result.alternatives = alternatives or [Alternative('A', -100.0, dict(amounts), 3, ['initial'])]
    return result


def two_result_sweep():
    b = {('gas', 'methane'): 0.1, ('gas', 'carbon dioxide'): 0.4, ('solid', 'carbon'): 0.5}
    second = fake_result(400.0)
    second.alternatives = [Alternative('A', -100.0, dict(second.amounts), 2, ['initial', 'random 0']),
                           Alternative('B', -99.5, b, 1, ['max gas'],
                                       phase_gibbs={'gas': -60.0, 'solid': -39.5})]
    points = [(300.0, 240.0, fake_result(300.0)), (400.0, 240.0, second),
              (450.0, 240.0, fake_result(450.0, dropped=['aqueous carbon'])), (500.0, 240.0, 'failed here')]
    return SweepResult(case(temperature_C=[300, 400, 450, 500]), points)


def test_summary_has_the_new_columns():
    summary = two_result_sweep().summary()
    for column in ('n_results', 'carbon_bound_violation', 'notes'):
        assert column in summary
    assert list(summary.n_results[:3]) == [1, 2, 1]
    assert summary.notes[2] == 'aqueous carbon bound dropped' and summary.notes[0] == ''


def test_alternatives_table_and_choices(tmp_path):
    sweep = two_result_sweep()
    folder = sweep.save(tmp_path / 'run')
    alternatives = pd.read_csv(folder / 'alternatives.csv')
    assert list(alternatives.columns[:8]) == ['T_C', 'P_bar', 'result', 'gibbs_J_per_mol_C', 'n_starts',
                                              'phase', 'species', 'mol_per_mol_C']
    assert set(alternatives.T_C) == {400.0}  # only the point with more than one result
    assert sorted(set(alternatives.result)) == ['A', 'B']
    assert not (folder / 'choices.toml').exists()
    tables = load_saved(folder)
    assert tables.flagged_points == [(400.0, 240.0)] and tables.choices == {}
    row = tables.summary[tables.summary.T_C == 400.0].iloc[0]
    assert row.gibbs_J_per_mol_C == -100.0 and row.result == 'A'

    write_choices(folder, {(400, 240): 'B', (300, 240): 'A'})
    assert read_choices(folder) == {(400.0, 240.0): 'B'}
    tables = load_saved(folder)
    assert tables.choices == {(400.0, 240.0): 'B'}
    row = tables.summary[tables.summary.T_C == 400.0].iloc[0]
    assert row.gibbs_J_per_mol_C == -99.5 and row.result == 'B'
    assert row.gas_carbon_percent == pytest.approx(50.0)
    species = tables.species[tables.species.T_C == 400.0].set_index(['phase', 'species'])
    assert species.loc[('solid', 'carbon'), 'mol_per_mol_C'] == 0.5
    assert species.loc[('solid', 'carbon'), 'mol_per_kg_dry_feed'] == pytest.approx(0.5 * 30)
    phases = tables.phases[tables.phases.T_C == 400.0].set_index('phase')
    assert phases.loc['solid', 'carbon_percent'] == pytest.approx(50.0)
    assert phases.loc['solid', 'gibbs_J_per_mol_C'] == pytest.approx(-39.5)
    # other points are untouched; going back to A restores the saved tables
    assert tables.summary[tables.summary.T_C == 300.0].iloc[0].gibbs_J_per_mol_C == -100.0
    back = tables.with_choices({})
    pd.testing.assert_frame_equal(back.species.reset_index(drop=True),
                                  load_saved(folder).with_choices({}).species.reset_index(drop=True))
    assert back.summary[back.summary.T_C == 400.0].iloc[0].gibbs_J_per_mol_C == -100.0
    write_choices(folder, {(400.0, 240.0): 'A'})
    assert not (folder / 'choices.toml').exists()


def test_apply_choices_ignores_unknown_results():
    sweep = two_result_sweep()
    summary, species, phases = sweep.summary(), sweep.species_table(), sweep.phase_summary()
    same = apply_choices(summary, species, phases, sweep.alternatives_table(), {(300.0, 240.0): 'C'})
    pd.testing.assert_frame_equal(same[1], species)


def test_old_saved_results_still_load(tmp_path):
    folder = tmp_path / 'old'
    folder.mkdir()
    pd.DataFrame({'T_C': [300.0], 'P_bar': [240.0], 'status': ['solved'], 'gibbs_J_per_mol_C': [-1.0],
                  'gas_carbon_percent': [50.0], 'error': [np.nan]}).to_csv(folder / 'summary.csv', index=False)
    pd.DataFrame({'T_C': [300.0], 'P_bar': [240.0], 'phase': ['gas'], 'species': ['methane'],
                  'mol_per_mol_C': [0.5]}).to_csv(folder / 'species.csv', index=False)
    pd.DataFrame({'T_C': [300.0], 'phase': ['gas'], 'carbon_percent': [50.0]}).to_csv(folder / 'phases.csv',
                                                                                     index=False)
    (folder / 'case.toml').write_text((Path(__file__).parent / 'data' / 'legacy_cases' /
                                       'case_template.toml').read_text())
    tables = load_saved(folder)
    assert tables.alternatives.empty and tables.choices == {} and tables.flagged_points == []
    assert list(tables.summary.columns) == ['T_C', 'P_bar', 'status', 'gibbs_J_per_mol_C',
                                            'gas_carbon_percent', 'error']


def test_live_tables_expose_alternatives(tmp_path):
    sweep = two_result_sweep()
    run = SimpleNamespace(name='x', id=1, case=sweep.case, output_dir=tmp_path, sweep_result=lambda: sweep)
    write_choices(tmp_path, {(400.0, 240.0): 'B'})
    tables = live_tables(run)
    assert tables.flagged_points == [(400.0, 240.0)]
    assert tables.summary[tables.summary.T_C == 400.0].iloc[0].result == 'B'


# ---------------------------------------------------------------- solving (Julia)
def test_cutoff_case_solves_at_every_point(base_model):
    """aqueous = {exact = 0.03} with the cut-off inside the range: every point solves,
    and the points at and above the cut-off say that the bound was dropped."""
    c = case(temperature_C=[350, 400], carbon_bounds={'aqueous': {'exact': 0.03}}, aqueous_cutoff_C=390,
             continuation=False)
    sweep = base_model.sweep(c, progress=None)
    summary = sweep.summary()
    assert list(summary.status) == ['solved'] * 3, summary.error.tolist()
    assert list(summary.T_C) == [350.0, 390.0, 400.0]
    assert list(summary.notes) == ['', 'aqueous carbon bound dropped', 'aqueous carbon bound dropped']
    below = sweep.points[0][2]
    aqueous_carbon = sum(v * below.species_data[k]['elements'].get('C', 0) for k, v in below.amounts.items()
                         if k[0] == 'aqueous')
    assert aqueous_carbon / below.feed_atoms[below.elements.index('C')] == pytest.approx(0.03, abs=1e-8)


def test_general_bounds_are_respected(base_model):
    c = case(temperature_C=450, carbon_bounds={'gas': {'min': 0.5, 'max': 0.9}, 'solid': {'max': 0.02}},
             aqueous_cutoff_C=390)
    r = base_model.solve(c)
    feed_C = r.feed_atoms[r.elements.index('C')]
    share = {}
    for (phase, s), v in r.amounts.items():
        share[phase] = share.get(phase, 0.0) + v * r.species_data[(phase, s)]['elements'].get('C', 0) / feed_C
    assert 0.5 - 1e-7 <= share['gas'] <= 0.9 + 1e-7
    assert share['solid'] <= 0.02 + 1e-7
    assert r.carbon_bound_violation < 1e-7
