"""Species discovery must find exactly the species (and data) of the original model."""
import pandas as pd
import pytest

from helpers import same

KINDS = ['base', 'acetone']


def records(frame):
    return frame.astype(object).where(frame.notna(), None).to_dict('records')


@pytest.mark.parametrize('kind', KINDS)
def test_catalogue(models, catalogue_ref, kind):
    new = records(models[kind].species(include_unavailable=True))
    ref = catalogue_ref[kind]['catalogue']
    assert len(new) == len(ref)
    for a, b in zip(new, ref):
        assert same(a, b), (a, b)


@pytest.mark.parametrize('kind', KINDS)
@pytest.mark.parametrize('registry', ['molecular', 'organic', 'solid', 'ion'])
def test_registries(models, catalogue_ref, kind, registry):
    reg = getattr(models[kind].registry, registry)
    ref = catalogue_ref[kind]['registries'][registry]
    assert list(reg) == list(ref)
    for name, expected in ref.items():
        actual = {k: reg[name].get(k) for k in expected}
        assert same(actual, expected), name


@pytest.mark.parametrize('kind', KINDS)
@pytest.mark.parametrize('table', ['ion_epcsaft_table', 'ion_hkf_table'])
def test_ion_tables(models, catalogue_ref, kind, table):
    new = records(getattr(models[kind].registry, table))
    assert same(new, catalogue_ref[kind][table])


def test_mineral_inventory(base_model, catalogue_ref):
    assert same(records(base_model.mineral_inventory()), catalogue_ref['base']['mineral_inventory'])


def test_product_candidates(base_model, catalogue_ref):
    # The original default element set; source labels were renamed, so they are not compared.
    new = base_model.product_candidates(['Na', 'K', 'Cs', 'Li', 'Cl', 'Mg'], temperature_C=400.0)
    ref = pd.DataFrame(catalogue_ref['base']['ion_product_candidates_400C'])
    columns = ['Species', 'Phase', 'Elements', 'Ready', 'Reason']
    assert same(records(new[columns]), records(ref[columns]))
