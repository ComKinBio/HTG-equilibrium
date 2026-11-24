# HTG equilibrium

Chemical equilibrium for hydrothermal gasification (HTG) by constrained
Gibbs-energy minimization:

- fluid phases (gas, organic liquid, molecular aqueous) with PC-SAFT, via
  [Clapeyron.jl](https://github.com/ClapeyronThermo/Clapeyron.jl)
- an aqueous electrolyte phase with ePC-SAFT and SUPCRT92/HKF ion standard states
- pure solids from NASA-7 polynomials and SUPCRT92 minerals
- element, charge and optional carbon-allocation constraints
- IPOPT (default), SLSQP or trust-constr optimizers, with several starting
  points, detection of distinct local minima, and optimality diagnostics

## Setup

1. Install [uv](https://docs.astral.sh/uv/getting-started/installation/).
2. In the repository folder:

   ```sh
   uv sync
   ```

The first calculation downloads Julia (if it is not already installed) and
installs the Julia packages into `.venv/julia_env`. This takes several minutes
once, but only needs to be done once.

## The app

```sh
uv run htg gui
```

opens the app in its own window (`uv run htg gui --browser`: in a browser tab).
It runs on your own computer (only you can open it) until you close the window,
or the terminal, or press Ctrl+C. Closing the window stops the workers.

The window uses Qt (installed by `uv sync`). On Linux it also needs a desktop
session and a few system libraries; on Debian or Ubuntu:

```sh
sudo apt install libegl1 libxcb-cursor0 libxkbcommon-x11-0 libnss3
```

If the window cannot open, the app prints one line saying why and opens in the
browser instead. Tabs:

- **Case**: open, edit and save case files (comments in the file are kept).
  Temperatures and pressures can be entered in °C or K and bar, MPa, Pa or atm.
  The feed is a total mass, solids, ash and an elemental analysis of the
  ash-free solids (or amounts of species), with a live summary of the implied
  masses and moles. Carbon bounds per phase, solver settings, and *Advanced
  settings* follow. *Check* runs the fast checks; *Add to queue* starts the
  calculation.
- **Queue**: runs in progress, and the number of parallel workers. Workers
  start when the app opens; each needs about a minute and 1.5 GB of memory
  before its first point. All queued runs share the workers.
- **Results**: plots of amounts (all phases together or per phase, per kg of
  total, dry or dry ash-free feed, in mol or g), gas composition, carbon
  gasification efficiency, carbon in each phase and temperature–pressure maps,
  with tables of this session's runs and of saved results. *Save plot…* writes
  PNG, SVG or PDF at a chosen size. When different starting points ended at
  different minima, a banner offers to compare them and to choose which result
  the plots and tables use.
- **Compare**: one quantity from several runs in the same plot, and the
  settings in which the runs differ.
- **Species**: every species in the databases, filters by element, data source
  and phase, and buttons that add a species to the case being edited. *New
  species*, *Edit selected* and *Import from file…* write to your own species
  database (below).

Display choices (units, colours, line width, export sizes) are kept in
`settings.toml` in your config folder (`~/.config/htg-equilibrium/` on Linux).
They never change a calculation.

Every run is saved to `results/<case name>-<date>-<time>/`. Options:
`--browser`, `--workers N`, `--results DIR`, `--cases DIR`, `--port PORT`,
`--no-browser` (only start the server).

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
| `alternatives.csv` | points where different starting points ended at different minima: every result (A, B, …) with its Gibbs energy and amounts |
| `choices.toml` | written by the app: which result each point uses in plots and tables (missing: the lowest Gibbs energy, A) |
| `user_species_used.toml` | the records from your own species database that the run used (only when it used any) |
| `species.png`, `carbon.png` | species profiles and carbon allocation against temperature |
| `case.toml` | the case that produced the results |

`summary.csv` also has `n_results` (number of distinct results at the point),
`carbon_bound_violation` (largest violation of a carbon bound, as a fraction of
feed carbon) and `notes` (for example when the aqueous carbon bound was dropped
above the aqueous cut-off).

Temperatures and pressures can both be lists (every temperature is solved at
every pressure), or a case can list explicit `points = [[T, P], ...]`.

### Case files

Case files use format 2 (`format = 2`); `htg new` writes a commented template.
Old case files still load: the app converts them when you open them and
writes format 2 when you save. The main parts:

```toml
format = 2

[feed]
total_kg = 10.0                           # solids + water
solids = { value = 15.0, unit = "%" }     # "%" of the total, or "kg"
ash = { value = 8.45, unit = "%" }        # "%" of the solids, or "kg"; leave out for no ash
elements_wt_percent = { C = 30.4, H = 4.41, O = 53.58, N = 2.66, S = 0.51 }
# % of the ash-free solids; scaled to 100 % if they do not add up

[constraints.carbon]                      # share of feed carbon per phase (fractions)
gas = { min = 0.90 }                      # exact, or min and/or max
solid = { max = 0.05 }

[solver]
starts = 5                                # starting points per point, in total
```

Instead of the elemental analysis, `[feed.species]` gives amounts of species
(`water = { value = 4.5, unit = "kg" }`, units `kg` or `mol`). A carbon bound
on the aqueous phase is dropped at temperatures at and above the aqueous
cut-off; the point is solved and the summary notes it.

In the app, temperature and pressure fields accept lists and ranges:

| Format | Type | Gives |
|---|---|---|
| List | `300, 350, 400` | 300, 350, 400 |
| Range (start:stop:step) | `300:400:50` | 300, 350, 400 |
| Evenly spaced (start:stop n=count) | `300:400 n=5` | 300, 325, 350, 375, 400 |
| Mixed | `250, 300:400:50` | 250, 300, 350, 400 |

### Starting points and distinct results

The optimizer finds local minima of the Gibbs energy. With `starts = N` each
point is solved from N starting points. Two solutions count as different
results when any amount differs by more than `distinct_composition_tol`
(mol/mol C) or the Gibbs energies differ by more than `distinct_gibbs_tol`
(J/mol C). The results are labelled A, B, C… by increasing Gibbs energy; A is
used unless you choose another one in the app.

Options:

- `htg run a.toml b.toml` runs several cases.
- `--workers N` solves in N parallel processes. Independent pieces of all cases
  (single points, or continuation chains: one per pressure) share the workers.
  Each worker first starts Julia (about a minute, longer when many start at
  once), so this pays off for large grids.
- `--out DIR` chooses the output folder.
- `htg plot results/my_case` redraws the plots from the saved tables.
- `--user-data DIR` uses another folder for your own species database;
  `--no-user-species` uses the shipped databases only (both for `run`, `check`,
  `species` and `gui`).

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
| `cases/` | example case files (format 2) and comparison data |
| `docs/` | guides for extending the model |
| `tests/` | tests, including reference results of the original implementation |

Package modules, from inputs to results:

| Module | Role |
|---|---|
| `case.py` | `Case` and `Feed`: case definition, TOML case files, conversion of old files |
| `feed_v2.py` | the feed as entered (total, solids, ash, elements or species) and its summary |
| `databases.py` | reading the database files |
| `species.py` | species registry: `species.toml` plus discovery of every species with complete data |
| `user_species.py`, `species_import.py`, `species_catalogue.py` | your own species database, import from files, the catalogue and its filters |
| `saft.py` | looking up PC-SAFT/ePC-SAFT parameter sets |
| `standard_states.py` | standard-state Gibbs energies (NASA-7, SUPCRT minerals, HKF ions) |
| `phases.py` | phase models: Gibbs energy and chemical potentials of each phase |
| `equilibrium.py` | balances, constraints, feasible start, and the solve |
| `optimizer.py` | optimizer backends, starting points, grouping into distinct results, optimality diagnostics |
| `sweep.py` | temperature-pressure sweeps, preflight checks |
| `runner.py` | pool of worker processes with a shared queue of pieces |
| `results.py`, `plotting.py` | result tables, alternatives and choices, plots for `htg run` |
| `units.py`, `settings.py` | unit conversion and display; app preferences |
| `quantities.py` | amounts per kg of feed, gas composition, carbon gasification efficiency |
| `plotspec.py`, `plotbuild.py`, `plotrender.py` | plots as plain data, built from result tables, drawn with Plotly (app) or Matplotlib (exports) |
| `gui/` | the app: one module per tab (`editor.py` is the Case tab), dialogs, the window (`window.py`) and the page (`app.py`) |
| `model.py` | `Model`: databases + registry, the entry point |
| `cli.py` | the `htg` command |

## Extending the model

- [Adding species](docs/adding_species.md): in the app (Species tab), or by
  adding database rows.
- [Adding a phase model](docs/adding_phase_model.md): a new equation of state
  or activity model for a phase.

## Your own species database

Species you add or edit in the app (Species tab: *New species*, *Edit
selected*, *Import from file…*) are saved in your own species database, a
layer over `data/`:

- Linux: `~/.local/share/htg-equilibrium/species/` (other systems: the user
  data folder for `htg-equilibrium`)
- another folder: `--user-data DIR` or the environment variable `HTG_USER_DATA`

It holds `user_pcsaft.csv`, `user_epcsaft.csv` (same columns as the shipped
files), `user_thermo.txt` (NASA-7 records), `user_supcrt.json` (HKF and
Maier–Kelley records) and `user_species.toml`. A record with the same species
name as a shipped one replaces it. The shipped databases in `data/` are never
changed. Every saved run lists the user records it used in
`user_species_used.toml`, and its `case.toml` notes the folder.

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
