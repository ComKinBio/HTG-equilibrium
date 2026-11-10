"""Standard states, SAFT parameters and phase models at fixed states and compositions.

These compare the building blocks of the Gibbs energy directly, independently
of the optimizer.
"""
import numpy as np
import pytest

from htg import standard_states as ss
from htg.constants import ATOMIC_WEIGHTS
from htg.databases import nasa7_properties
from htg.phases import ElectrolytePhase, FluidPhase, build_pcsaft_model

from helpers import close, same
from convert import load

STATES = load('thermo')['states']
IDS = [f"{s['T_C']:g}C-{s['P_bar']:g}bar" for s in STATES]
PARAM_KEYS = ('Formula', 'Source', 'Database row', 'Model', 'RefID', 'm', 'sigma', 'epsilon',
              'epsilonAB', 'kappaAB', 'association_scheme', 'Mw')


def kelvin_pascal(state):
    return state['T_C'] + 273.15, state['P_bar'] * 1e5


def outcome(function):
    try:
        return function()
    except Exception as exc:
        return f'{type(exc).__name__}: {exc}'


@pytest.mark.parametrize('state', STATES, ids=IDS)
def test_molecular_g0(base_model, state):
    T, _ = kelvin_pascal(state)
    reg, db = base_model.registry, base_model.databases.nasa
    assert sorted(reg.all_molecular) == list(state['molecular_g0'])
    for name, expected in state['molecular_g0'].items():
        actual = outcome(lambda: nasa7_properties(db, reg.nasa_name(name), T)[3])
        assert same(actual, expected), name


@pytest.mark.parametrize('state', STATES, ids=IDS)
def test_solid_g0(base_model, state):
    T, P = kelvin_pascal(state)
    reg, db = base_model.registry, base_model.databases.nasa
    assert sorted(reg.solid) == list(state['solid_g0'])
    for name, expected in state['solid_g0'].items():
        actual = outcome(lambda: float(ss.solid_g0(reg, db, [name], T, P)[0]))
        assert same(actual, expected), name


@pytest.mark.parametrize('state', STATES, ids=IDS)
def test_saft_parameters(base_model, state):
    T, _ = kelvin_pascal(state)
    for name, expected in state['saft_parameters'].items():
        def get():
            d = base_model.registry.resolve_parameters(base_model.databases, T, [name])[name]
            return {k: d.get(k) for k in PARAM_KEYS}
        assert same(outcome(get), expected), name


FLUID_CASES = [(s, label) for s in STATES for label in s['fluid_phases']]


@pytest.mark.parametrize('state,label', FLUID_CASES,
                         ids=[f"{s['T_C']:g}C-{label}" for s, label in FLUID_CASES])
def test_fluid_phase(base_model, state, label):
    ref = state['fluid_phases'][label]
    T, P = kelvin_pascal(state)
    reg, names = base_model.registry, ref['species']
    mu0 = ss.nasa_g0(reg, base_model.databases.nasa, names, T)
    close(mu0, ref['mu0'])
    resolved = reg.resolve_parameters(base_model.databases, T, names)
    phase = FluidPhase(label, names, mu0, build_pcsaft_model(names, resolved, reg.binary_interactions),
                       T, P, ref['root'])
    for point in ref['points']:
        n = np.array([point['amounts'][s] for s in names])
        if 'error' in point:
            with pytest.raises(ValueError):
                phase.gibbs(n)
            continue
        assert phase.gibbs(n) == pytest.approx(point['G'], rel=1e-10, abs=1e-8)
        close(phase.lnphi(n), point['ln_phi'], rtol=1e-9, atol=1e-11)


ELECTROLYTE_CASES = [(s, label) for s in STATES for label in s.get('electrolyte', {})]


@pytest.mark.parametrize('state,label', ELECTROLYTE_CASES,
                         ids=[f"{s['T_C']:g}C-{label}" for s, label in ELECTROLYTE_CASES])
def test_electrolyte_phase(base_model, state, label):
    ref = state['electrolyte'][label]
    T, P = kelvin_pascal(state)
    reg = base_model.registry
    selection = ref['selection']
    neutrals = list(dict.fromkeys(reg.aqueous_reference_neutrals
                                  + [s for s in selection if s not in reg.ion]))
    resolved = reg.resolve_parameters(base_model.databases, T, neutrals, ATOMIC_WEIGHTS)
    phase = ElectrolytePhase(base_model, selection, T, P, resolved, ATOMIC_WEIGHTS)
    assert phase.all_species == ref['aqueous_species']
    assert same(phase.water, ref['water'])
    close(phase.mu0_ion, ref['mu0_ion'])
    close(phase.mu0_neutral, ref['mu0_neutral'])
    close(phase.lnphi_infinity, ref['lnphi_infinity'], rtol=1e-9, atol=1e-11)
    for point in ref['points']:
        full = np.array([point['amounts'].get(s, 0.0) for s in phase.all_species])
        if 'error' in point:
            with pytest.raises(ValueError):
                phase.gibbs_full(full)
            continue
        assert phase.gibbs_full(full) == pytest.approx(point['G'], rel=1e-10, abs=1e-8)
        close(np.log(phase.phi(full)), point['ln_phi'], rtol=1e-9, atol=1e-11)
