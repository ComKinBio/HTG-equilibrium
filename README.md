# HTG equilibrium

Chemical equilibrium for hydrothermal gasification (HTG) by constrained
Gibbs-energy minimization:

- fluid phases (gas, organic liquid, molecular aqueous) with PC-SAFT, via
  [Clapeyron.jl](https://github.com/ClapeyronThermo/Clapeyron.jl)
- an aqueous electrolyte phase with ePC-SAFT and SUPCRT92/HKF ion standard states
- pure solids from NASA-7 polynomials and SUPCRT92 minerals
- element, charge and optional carbon-allocation constraints
- IPOPT (default), SLSQP or trust-constr optimizers, with multistart and
  optimality diagnostics

## Setup

1. Install [uv](https://docs.astral.sh/uv/getting-started/installation/).
2. In the repository folder:

   ```sh
   uv sync
   ```

The first calculation downloads Julia (if it is not already installed) and
installs the Julia packages into `.venv/julia_env`. This takes several minutes
once; later runs start in seconds.

## The app

```sh
uv run htg gui
```

opens the app in your web browser. It runs on your own computer (only you can
open it) until you close the terminal or press Ctrl+C. Tabs:

- **Case**: open, edit and save case files (comments in the file are kept).
  Basic settings are shown first; the rest is under *Advanced settings*.
  *Check* runs the fast checks; *Add to queue* starts the calculation.
- **Queue**: runs in progress, and the number of parallel workers. Workers
  start when the app opens; each needs about a minute and 1.5 GB of memory
  before its first point. All queued runs share the workers.
- **Results**: interactive plots (species against temperature or pressure,
  temperature–pressure maps, carbon allocation) and tables of this session's
  runs and of saved results. Use the camera icon on a plot to save a PNG.
- **Compare**: one species or phase across several runs, and the settings in
  which the runs differ.
- **Species**: every species in the databases, why unavailable ones cannot be
  used, and buttons that add a species to the case being edited.

Every run is saved to `results/<case name>-<date>-<time>/`. Options:
`--workers N`, `--results DIR`, `--cases DIR`, `--port PORT`, `--no-browser`.

## Running a case from the terminal

A *case file* describes one calculation or a temperature sweep: conditions,
feed, allowed products in each phase, constraints and solver. Examples are in
[`cases/`](cases/).

```sh
uv run htg new my_case.toml          # commented template to edit
uv run htg species --search ethanol  # names of species you can use
uv run htg check my_case.toml        # fast checks, no equation of state
uv run htg run my_case.toml          # solve; writes results/my_case/
```

`htg run` writes to `results/<case name>/`:

| File | Content |
|---|---|
| `summary.csv` | one row per point: status, Gibbs energy, gas-carbon share, errors |
| `species.csv` | amount of every species: mol/mol feed C, mol and g per kg dry feed, molality of ions |
| `phases.csv` | per phase: amount, mass, share of feed carbon, Gibbs energy |
| `species.png`, `carbon.png` | species profiles and carbon allocation against temperature |
| `case.toml` | the case that produced the results |

Temperatures and pressures can both be lists (every temperature is solved at
every pressure), or a case can list explicit `points = [[T, P], ...]`.

Options:

- `htg run a.toml b.toml` runs several cases.
- `--workers N` solves in N parallel processes. Independent pieces of all cases
  (single points, or continuation chains: one per pressure) share the workers.
  Each worker first starts Julia (about a minute, longer when many start at
  once), so this pays off for large grids.
- `--out DIR` chooses the output folder.
- `htg plot results/my_case` redraws the plots from the saved tables.

The example cases reproduce the original notebooks:

| Case | Source |
|---|---|
| `yakaboylu2013_fig7_8.toml` | Swine manure, Yakaboylu et al. (2013) Figures 7–8 |
| `cellulose_matlab.toml` | Cellulose compared with the former MATLAB model (`data/Cellulose.xlsx`) |
| `manure_ions_minerals.toml` | Manure with aqueous ions and minerals, 80 % gas carbon |
| `manure_gas_min.toml` | As above with salts and at least 90 % gas carbon |
| `manure_configurable.toml` | Default case of the former configurable notebook |

## Using Python

```python
import htg

model = htg.Model.load()                       # reads data/
case = htg.load_case('cases/yakaboylu2013_fig7_8.toml')

result = model.solve(case, temperature_C=450)  # one temperature
print(result)
result.species_table()                         # pandas DataFrame

sweep = model.sweep(case)                      # all temperatures of the case
sweep.summary()
sweep.save('results/fig78')

case = case.replace(pressure_bar=300, solver='slsqp')   # modified copy
```

## Repository layout

| Path | Content |
|---|---|
| `src/htg/` | the package (see below) |
| `data/` | databases: NASA-7 (`Thermdat_WasteMech.txt`), PC-SAFT and ePC-SAFT parameters (`PROJECT_*.csv`), SUPCRT92 (`supcrt92.xml`), species registry (`species.toml`) |
| `cases/` | example case files and comparison data |
| `docs/` | guides for extending the model |
| `tests/` | tests, including reference results of the original implementation |

Package modules, from inputs to results:

| Module | Role |
|---|---|
| `case.py` | `Case` and `Feed`: case definition, TOML case files |
| `databases.py` | reading the database files |
| `species.py` | species registry: `species.toml` plus discovery of every species with complete data |
| `saft.py` | looking up PC-SAFT/ePC-SAFT parameter sets |
| `standard_states.py` | standard-state Gibbs energies (NASA-7, SUPCRT minerals, HKF ions) |
| `phases.py` | phase models: Gibbs energy and chemical potentials of each phase |
| `equilibrium.py` | balances, constraints, feasible start, and the solve |
| `optimizer.py` | optimizer backends, multistart, optimality diagnostics |
| `sweep.py` | temperature-pressure sweeps, preflight checks |
| `runner.py` | pool of worker processes with a shared queue of pieces |
| `results.py`, `plotting.py` | result tables and plots |
| `gui/` | the app: case editor (`editor.py`, `document.py`), result plots (`plots.py`), pages (`app.py`) |
| `model.py` | `Model`: databases + registry, the entry point |
| `cli.py` | the `htg` command |

## Extending the model

- [Adding species](docs/adding_species.md): molecules, ions and solids usually
  need only database rows.
- [Adding a phase model](docs/adding_phase_model.md): a new equation of state
  or activity model for a phase.

## Tests

```sh
uv run pytest               # about 4 minutes
uv run pytest -m "not slow"
```

The thermodynamic model was restructured from notebooks into this package. The
tests check that the package reproduces the original implementation:
standard-state values, SAFT parameters, phase Gibbs energies and fugacity
coefficients at fixed compositions, and complete equilibrium solves. The
reference values were generated by a frozen copy of the original code; see
[`tests/reference/README.md`](tests/reference/README.md). Any change that alters
the thermodynamics makes these tests fail.
