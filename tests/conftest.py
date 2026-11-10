"""Shared fixtures. The reference data in tests/reference/data come from the frozen
original implementation (see tests/reference/README.md)."""
import sys
from pathlib import Path

import pytest

REFERENCE = Path(__file__).resolve().parent / 'reference'
# The databases and species.toml the reference results were generated with. The
# reference tests use this frozen copy, so editing data/ (adding species) never
# breaks them; they test that the code reproduces the original model.
FROZEN_DATA = REFERENCE / 'databases'
sys.path.insert(0, str(REFERENCE))

import htg  # noqa: E402
from convert import load  # noqa: E402


@pytest.fixture(scope='session')
def models():
    """'base': registry without extra NASA names; 'acetone': the shipped species.toml."""
    return {'base': htg.Model.load(FROZEN_DATA, nasa_names={}),
            'acetone': htg.Model.load(FROZEN_DATA)}


@pytest.fixture(scope='session')
def base_model(models):
    return models['base']


@pytest.fixture(scope='session')
def catalogue_ref():
    return load('catalogue')


@pytest.fixture(scope='session')
def thermo_ref():
    return load('thermo')


@pytest.fixture(scope='session')
def solves_ref():
    return load('solves')
