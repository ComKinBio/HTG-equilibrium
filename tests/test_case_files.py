"""The case files in cases/ encode the same problems as the original notebooks.

Each case file, at a reference temperature, must give the same balances,
feasible start, bounds and Gibbs energy function as the corresponding
reference case (built from the notebook settings).
"""
from pathlib import Path

import numpy as np
import pytest

import htg
from htg.equilibrium import setup

from convert import case_from_legacy, load, pairs

CASES = Path(__file__).resolve().parents[1] / 'cases'
REFERENCE = {c['name']: c for c in load('solves')['cases']}
MATCHES = [
    ('yakaboylu2013_fig7_8', 'fig78_300C'), ('yakaboylu2013_fig7_8', 'fig78_380C'),
    ('yakaboylu2013_fig7_8', 'fig78_450C'),
    ('cellulose_matlab', 'cellulose_400C'), ('cellulose_matlab', 'cellulose_500C'),
    ('manure_ions_minerals', 'full_text_250C'), ('manure_ions_minerals', 'full_text_400C'),
    ('manure_gas_min', 'full_text_gas_min_300C'),
    ('manure_configurable', 'configurable_300C'), ('manure_configurable', 'configurable_450C'),
]


def test_all_case_files_load():
    files = sorted(CASES.glob('*.toml'))
    assert files
    for path in files:
        case = htg.load_case(path)
        assert case.name == path.stem
        if case.comparison:
            assert Path(case.comparison['file']).exists()


def test_round_trip(tmp_path):
    for path in sorted(CASES.glob('*.toml')):
        case = htg.load_case(path)
        case.save(tmp_path / path.name)
        again = htg.load_case(tmp_path / path.name)
        again.comparison = case.comparison  # paths are resolved relative to the file
        assert again == case, path.name


def test_template_loads():
    template = Path(htg.__file__).parent / 'templates' / 'case_template.toml'
    htg.load_case(template)


@pytest.mark.parametrize('file,reference', MATCHES, ids=[f'{f}-{r}' for f, r in MATCHES])
def test_case_file_matches_reference(models, file, reference):
    entry = REFERENCE[reference]
    model = models['acetone']  # the shipped species.toml
    T = entry['settings']['temperature_C']
    new = setup(model, htg.load_case(CASES / f'{file}.toml').at(T))
    old = setup(model, case_from_legacy(entry))
    assert new.pairs == old.pairs
    assert new.elements == old.elements
    for attr in ('b', 'A_full', 'b_full', 'initial_amounts', 'x_start'):
        np.testing.assert_allclose(getattr(new, attr), getattr(old, attr), rtol=1e-12, atol=1e-14,
                                   err_msg=attr)
    assert new.bounds == old.bounds
    assert new.gas_min_fraction == old.gas_min_fraction
    final = np.array(list(pairs(entry['result']['amounts']).values()))
    assert new.total_G(final) == pytest.approx(old.total_G(final), rel=1e-13)
