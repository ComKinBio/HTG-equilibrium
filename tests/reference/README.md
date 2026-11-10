# Reference results

The thermodynamic model in `src/htg` must reproduce the original notebook
implementation. This folder pins that implementation:

- `legacy/` holds a **frozen** copy of the model code (`HTG_Model.ipynb` code
  cells) and `htg_optimizer.py` as of commit `74ac11d`. Never edit it.
- `legacy_cases.py` defines reference cases (from the original case notebooks)
  in the legacy settings format, plus the states and compositions used for
  composition-level checks.
- `databases/` is a frozen copy of the databases and `species.toml` the
  references were generated with. The tests load the model from here, so
  editing `data/` (for example adding species) does not affect them.
- `generate.py` runs the frozen code and writes `data/*.json`:
  - `catalogue.json`: species discovered from the databases, and registries
  - `thermo.json`: standard chemical potentials, SAFT parameters, phase Gibbs
    energies and fugacity coefficients at fixed compositions, electrolyte
    reference terms, at several temperatures and pressures
  - `solves.json`: complete equilibrium solves, including expected failures

`tests/test_reference_*.py` compare the package against these files. To add a
reference case, add it to `legacy_cases.py` and rerun
`uv run python tests/reference/generate.py` (about 3 minutes).
