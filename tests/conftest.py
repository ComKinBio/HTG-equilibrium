"""Shared fixtures. The reference data in tests/reference/data come from the frozen
original implementation (see tests/reference/README.md)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent / 'reference'))

import htg  # noqa: E402
from convert import load  # noqa: E402


@pytest.fixture(scope='session')
def models():
    """'base': registry without extra NASA names; 'acetone': the shipped species.toml."""
    return {'base': htg.Model.load(nasa_names={}), 'acetone': htg.Model.load()}


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
