"""Consistency of the phase models: mu = dG/dn and G homogeneous of degree one.

The optimizers use mu as the exact gradient of G; a new phase model should pass
the same checks (see docs/adding_phase_model.md).
"""
import numpy as np
import pytest

from htg import standard_states as ss
from htg.constants import ATOMIC_WEIGHTS
from htg.phases import ElectrolytePhase, FluidPhase, build_pcsaft_model

GAS = ['water', 'hydrogen', 'carbon dioxide', 'methane', 'nitrogen']
AQUEOUS = ['water', 'carbon dioxide', 'methane', 'Na+', 'HCO3-']
T, P = 300 + 273.15, 240e5


@pytest.fixture(scope='module')
def phases(base_model):
    reg, db = base_model.registry, base_model.databases
    resolved = reg.resolve_parameters(db, T, GAS)
    gas = FluidPhase('gas', GAS, ss.nasa_g0(reg, db.nasa, GAS, T),
                     build_pcsaft_model(GAS, resolved, reg.binary_interactions), T, P, 'vapour')
    neutrals = list(dict.fromkeys(reg.aqueous_reference_neutrals))
    resolved = reg.resolve_parameters(db, T, neutrals, ATOMIC_WEIGHTS)
    aqueous = ElectrolytePhase(base_model, AQUEOUS, T, P, resolved, ATOMIC_WEIGHTS)
    return {'gas': (gas, np.array([0.3, 0.05, 0.4, 0.2, 0.05])),
            'aqueous': (aqueous, np.array([10.0, 0.01, 0.002, 0.003, 0.003]))}


@pytest.mark.parametrize('name', ['gas', 'aqueous'])
def test_mu_is_gradient_of_G(phases, name):
    phase, n = phases[name]
    G, mu = phase.gibbs_and_mu(n)
    assert G == pytest.approx(phase.gibbs(n), rel=1e-12)
    for i, species in enumerate(phase.species):
        if species in ('Na+', 'HCO3-'):
            continue  # a single ion cannot be varied without breaking electroneutrality
        # Step large enough that round-off in G (~1e6 J) stays below the tolerance.
        h = 1e-4 * n[i]
        up, down = n.copy(), n.copy()
        up[i] += h
        down[i] -= h
        derivative = (phase.gibbs(up) - phase.gibbs(down)) / (2 * h)
        assert derivative == pytest.approx(mu[i], rel=1e-6, abs=1e-2), species


@pytest.mark.parametrize('name', ['gas', 'aqueous'])
def test_homogeneous(phases, name):
    phase, n = phases[name]
    assert phase.gibbs(2.5 * n) == pytest.approx(2.5 * phase.gibbs(n), rel=1e-10)
