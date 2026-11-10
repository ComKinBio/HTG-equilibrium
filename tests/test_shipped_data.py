"""The shipped data/ folder and case files work (without comparing to references)."""
from pathlib import Path

import pytest

import htg

CASES = sorted((Path(__file__).resolve().parents[1] / 'cases').glob('*.toml'))


@pytest.fixture(scope='module')
def shipped():
    return htg.Model.load()


@pytest.mark.parametrize('path', CASES, ids=[p.stem for p in CASES])
def test_case_file_preflight(shipped, path):
    table = shipped.check(htg.load_case(path))
    # The electrolyte model is invalid just below some cut-offs (low water density);
    # those points are expected to be flagged, everything else must be ready.
    bad = table[(table.status != 'ready') & ~table.error.str.contains('aqueous', case=False)]
    assert bad.empty, bad.to_string()
