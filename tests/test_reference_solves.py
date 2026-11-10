"""Complete equilibrium problems from the original case notebooks.

Two levels:
- setup: the same balances, feasible start and Gibbs energy function as the
  original (G evaluated at the original's solution must agree to ~1e-10);
- solve: the optimizer reaches the original's solution.
"""
import numpy as np
import pytest

from htg import SolveError
from htg.equilibrium import setup

from convert import LEGACY_PHASE, case_from_legacy, load, pairs

CASES = load('solves')['cases']
SOLVED = [c for c in CASES if 'result' in c]
FAILED = [c for c in CASES if 'error' in c]
SLOW_SECONDS = 10.0
REFERENCE = {c['name']: c for c in CASES}
# The original IPOPT settings stop early in these cases (the KKT residual of
# major species stays at several RT). The package then polishes the point with
# SLSQP and must reach the solution of the given (lower-G) reference case.
IMPROVED = {'full_text_250C_ipopt': 'full_text_250C'}


def params(cases):
    return [pytest.param(c, id=c['name'],
                         marks=[pytest.mark.slow] if c['seconds'] > SLOW_SECONDS else [])
            for c in cases]


def legacy_terms(terms):
    return {LEGACY_PHASE[p]: v for p, v in terms.items()}


@pytest.mark.parametrize('entry', params(SOLVED))
def test_setup_and_gibbs_function(models, entry):
    ref = entry['result']
    s = setup(models[entry['model']], case_from_legacy(entry))
    expected_pairs = list(pairs(ref['amounts']))
    assert s.pairs == expected_pairs
    np.testing.assert_allclose(s.initial_amounts, list(pairs(ref['initial_amounts']).values()),
                               rtol=1e-9, atol=1e-11)
    assert s.elements == ref['elements']
    np.testing.assert_allclose(s.b, ref['feed_atoms'], rtol=1e-12, atol=1e-14)
    assert s.C_mol_per_g == pytest.approx(ref['feed_C_mol_per_g_dry'], rel=1e-12)
    final = np.array(list(pairs(ref['amounts']).values()))
    terms = legacy_terms(s.phase_G(final))
    for phase, expected in ref['phase_gibbs_J_per_mol_feed_C'].items():
        assert terms[phase] == pytest.approx(expected, rel=1e-10, abs=1e-6), phase


@pytest.mark.parametrize('entry', params(SOLVED))
def test_solve(models, entry):
    r = models[entry['model']].solve(case_from_legacy(entry))
    if entry['name'] in IMPROVED:
        better = REFERENCE[IMPROVED[entry['name']]]['result']
        assert r.optimizer['polish']['used']
        assert r.gibbs_J_per_mol_C < entry['result']['gibbs_J_per_mol_feed_C']
        assert r.gibbs_J_per_mol_C == pytest.approx(better['gibbs_J_per_mol_feed_C'], rel=1e-8)
        # G is flat in the gas/aqueous water split here, so amounts agree less tightly than G.
        np.testing.assert_allclose(list(r.amounts.values()), list(pairs(better['amounts']).values()),
                                   rtol=1e-2, atol=1e-4)
        return
    if r.optimizer.get('method') == 'ipopt':
        assert not r.optimizer['polish']['used']
    ref = entry['result']
    assert r.gibbs_J_per_mol_C == pytest.approx(ref['gibbs_J_per_mol_feed_C'], rel=1e-9)
    expected = pairs(ref['amounts'])
    assert list(r.amounts) == list(expected)
    np.testing.assert_allclose(list(r.amounts.values()), list(expected.values()),
                               rtol=1e-5, atol=1e-7)
    for phase, value in ref['phase_gibbs_J_per_mol_feed_C'].items():
        assert legacy_terms(r.phase_gibbs_J_per_mol_C)[phase] == pytest.approx(value, rel=1e-8, abs=1e-3)
    assert r.gas_carbon_fraction == pytest.approx(ref['gas_carbon_fraction_actual'], abs=1e-7)
    assert r.balance_residual < 1e-7
    rules = ref['phase_allocation']
    assert r.phase_rules['aqueous_disabled'] == rules['aqueous_disabled']
    assert r.phase_rules['removed_aqueous_species'] == rules['removed_aqueous_species']
    assert r.phase_rules['initial_composition_reset'] == rules['initial_composition_reset']
    if 'optimizer' in ref and ref['optimizer'].get('method') != 'slsqp-fd':
        assert r.optimizer['method'] == ref['optimizer']['method']
        for new, old in zip(r.optimizer['starts'], ref['optimizer']['starts'], strict=True):
            assert new['start'].replace('_', ' ') == old['start']  # phase labels renamed
            assert (new['message'], new['nit']) == (old['message'], old['nit'])
            assert new['G'] == pytest.approx(old['G'], rel=1e-9)


@pytest.mark.parametrize('entry', params(FAILED))
def test_expected_failure(models, entry):
    with pytest.raises(SolveError) as info:
        models[entry['model']].solve(case_from_legacy(entry))
    # Same stage and cause; the advice at the end of some messages was reworded.
    expected = entry['error'].split(': ', 1)[1]
    first_sentence = expected.split('. ')[0]
    assert str(info.value).startswith(first_sentence), (str(info.value), expected)
