# HTGEquilibrium.jl

The Julia science core of the HTG equilibrium app (Phase 2, milestone M2 of
`docs/spec/julia_migration_plan.md`). It is ported from the Python package
`src/htg` one module at a time, keeping the Python names so the two are easy to
compare. So far: constants, formula parsing, water properties, HKF standard
states of aqueous ions, the SUPCRT92 failure-band check that the Julia GUI
will use to warn about sweeps, the data layer (databases, species
registry, user species folder, standard states, model loading), reading
case files (units, feed arithmetic), the phase models (G, mu, ln phi and
Hessians of the solid, fluid and electrolyte phases), the solver (setup of
one equilibrium point, IPOPT/SLSQP with multistart, the result data), the
result tables and CSV files, sweeps, the preflight check, the worker pool
on threads, and the command line (`htg run`, `check`, `species`, `new`; see
"Command line").

| File | Contents | Python counterpart |
|---|---|---|
| `src/constants.jl` | R, P0, reference states, atomic weights, phase names | `constants.py` |
| `src/formulas.jl` | `parse_formula`, `split_charge`, `molar_mass`, `formula_summary`, `hill_formula`, `parse_element_list` | `formulas.py` |
| `src/water.jl` | Water models: `SUPCRT92Water()` (original Fortran) and `IAPWS95Water()` (Clapeyron). `water_properties(model, T, P)` returns the same 8 quantities as `supcrt_water`. | `standard_states.supcrt_water` |
| `src/hkf.jl` | `HKFRecord`, g-function, revised HKF G, `hkf_ion_g0(model, records, T, P)` | `standard_states.hkf_ion_g0` |
| `src/supcrt92_bands.jl` | `supcrt92_failed`, `supcrt92_failure_bands`, `supcrt92_band_check` | – (new) |
| `vendor/supcrt92/H2O92D.f` | SUPCRT92 water code, unmodified | compiled inside pychnosz |
| `deps/build.jl` | Compiles `H2O92D.f` to `deps/build/libh2o92.so` | – |
| `test/runtests.jl` | Runs one test file per area | – |
| `test/test_formulas.jl` | Formula functions against Python outputs (`test/data/formulas_python.json`, written by `test/dump_formulas.py`) | – |
| `test/test_water.jl` | Golden tests, band check, IAPWS-95 self-consistency | `tests/test_reference_thermo.py` (electrolyte part) |
| `test/data/hkf_ions.toml` | HKF records of the 8 shipped ions, as the app's registry reads them | `registry.ion_hkf_table` |
| `test/dump_pychnosz.py` | Writes the fixture above and, with `--grid`, a wide pychnosz grid (not committed) | – |
| `src/pycompat.jl` | Python behaviour the port relies on: `ValueError` (abstract; the case/feed/unit errors subtype it), `py_repr`/`py_str`, `py_float`/`py_tryfloat`/`py_int` (with Python's messages), `py_round`, `py_format_g` (C/Python `%g`), `py_isclose`, `py_truth`, `py_sorted`, `isna`, `casefold`, `py_slice`, `read_text`, `py_splitlines` | – |
| `src/toml_ordered.jl` | `parse_toml`, `read_toml`: TOML as nested `OrderedDict`s in file order (as `tomllib`); used for case files, `species.toml` and the user species folder | `tomllib` |
| `src/units.jl` | Submodule `Units`: `to_file`, `from_file`, `format_number`, `display`, `field_label`, ... | `units.py` |
| `src/case.jl` | `Salt`, `Feed`, `Case`, `load_case`, `case_from_dict`, `validate`, `water_model_of`, `temperatures`, `pressures`, `grid`, `chains`, `at`, `replace`, carbon bounds, legacy keys, `legacy_to_v2`; writing: `to_toml`, `dump_case`, `save` | `case.py` |
| `src/feed_v2.jl` | `feed_v2_to_model`, `feed_summary`, `feed_fractions`, `elemental_masses`, `ash_percent`, `molar_mass_function`, `species_elements_function`, `registry_elements`, `legacy_feed_to_v2` | `feed_v2.py` |
| `src/databases.jl` | `Databases`, `load_databases`, NASA-7 parser (`read_nasa_lines`, `parse_nasa_header_elements`, `nasa7_properties`), SAFT CSV tables (`DataTable`), `supcrt92.xml` (`XMLNode`), `normalize_name` | `databases.py` |
| `src/saft.jl` | Parameter-set selection: `first_complete_epcsaft_set`, `first_complete_pcsaft_set`, `first_complete_ion_epcsaft_set`, `first_complete_hkf_set`, `normalize_ion_name`, `water_sigma` | `saft.py` |
| `src/registry.jl` | `Registry` (species discovery `discover!`, `all_molecular`, `solid`, `ion`, `ion_hkf_table`, `resolve_parameters`, catalogue of `CatalogueRow`s), `NO_SAFT_REASON`, `formula_elements`, `molecular_weight` | `species.py` |
| `src/user_species.jl` | Loading a user species folder: `resolve_user_dir`, `load_user_database`, `apply_user_layer!`, `hide_species!` | `user_species.py` (loading part) |
| `src/standard_states.jl` | `nasa_g0`, `supcrt_mineral_properties`, `solid_g0`, `hkf_records`, `hkf_ion_g0(water_model, hkf_table, ions, T, P)` | `standard_states.py` |
| `src/model.jl` | `Model`, `load_model` (databases + registry + user layer), `species`, `mineral_inventory`, `product_candidates` (gas rule: `nasa_phase_is`); `hkf_ion_g0(water_model, model, ...)`, `nasa_g0`, `solid_g0` of a `Model` | `model.py` (data part) |
| `test/test_units.jl`, `test/test_case.jl`, `test/test_feed.jl` | Ported Python tests; with `test/data/cases_python.json`, every case file and 119 variants against Python | `test_units.py`, `test_case_files.py`, `test_case_format.py` |
| `test/dump_cases.py`, `test/python_cases.jl` | Write / compare `test/data/cases_python.json` (~0.6 MB, not committed) | – |
| `test/data_setup.jl` | Shared setup of the three tests below: the reference models, Python's `same`/`close` | `tests/conftest.py`, `tests/helpers.py` |
| `test/test_databases.jl` | Readers and Python helpers (CSV dtypes, XML, NASA-7, ordered TOML) | – |
| `test/test_registry.jl` | `catalogue.json` (all of it), a user species folder, shipped `data/` smoke test, threaded reads | `tests/test_reference_catalogue.py`, `tests/test_user_species.py` (loading) |
| `test/test_standard_states.jl` | `thermo.json`: molecular and solid g°, SAFT parameters, mu0 of the fluid and electrolyte phases | `tests/test_reference_thermo.py` |
| `src/phases.jl` | `EOSFailure`, `PhaseError`, `association_sites`, `build_pcsaft_model`, `SolidPhase`, `FluidPhase`, `ElectrolytePhase` (water model selectable), `gibbs`, `gibbs_and_mu`, `mu_at`, `lnphi`, `hessian` (exact, ForwardDiff), `phi`, `gibbs_full`, `expand` | `phases.py` |
| `test/test_phases.jl` | `thermo.json`: fluid G and ln phi, electrolyte `lnphi_infinity`, G and ln phi; mu = dG/dn, homogeneity, Hessians against central differences, water found by name, water model, errors, threads | `tests/test_reference_thermo.py` (phase part), `tests/test_phases.py` |
| `test/data/user_species/`, `test/data/user_species_python.json` | A user species folder written by the Python app and what Python loads from it; written by `test/dump_user_species.py` | – |
| `src/results.jl` | `Alternative`, `Result` (same fields; plus `water_model`), `n_results`, `notes`, `for_alternative` | `results.py` (data part) |
| `src/optimizer.jl` | `GibbsProblem` (`expand`, `gibbs_and_mu`, `fun_and_grad`, `hessian_blocks`, `check_against_model`, `kkt`), `linprog` (HiGHS), `_ipopt` (Ipopt.jl), `_slsqp` (NLopt), `starting_points`, `candidate_starts`, `result_label`, `group_solutions`, `tangent_plane_test`, `_polish`, `optimize` | `optimizer.py` |
| `src/equilibrium.jl` | `SolveError`, `species_to_wt_percent`, `apply_phase_rules`, `Stage`, `Setup`, `setup`, `carbon_bound_rows`, `carbon_bound_violation`, `_optimize`, `solve` (also `Model.solve`) | `equilibrium.py`, `Model.solve` |
| `test/test_optimizer.jl` | Carbon rows, cut-off rule, grouping, starts, multistart with a stand-in problem, LP, backends, errors | `tests/test_backend.py` (optimizer part) |
| `test/test_equilibrium.jl`, `test/solves_setup.jl` | `solves.json` (setup, solves, failures, polish safeguard), solving checks, Result data, threads; `case_from_legacy` | `tests/test_reference_solves.py`, `tests/reference/convert.py`, `tests/test_backend.py` (solving part) |
| `test/dump_solves.py`, `test/data/solves_python_ipopt.json` | The Python app's IPOPT solution of every `solves.json` point (32 kB, committed) | – |
| `src/tables.jl` | Result tables as `DataTable`s (pandas-like: `DataTable(rows)`, `concat`, `setcolumn!`, `insert_column!`, `select`, `take_rows`), `csv_text`/`to_csv` (pandas `to_csv(index=False)`), `read_csv` (pandas `read_csv`), `to_string` (pandas `to_string(index=False)`, the printed tables) | pandas |
| `src/results.jl` (tables) | `species_table`, `phase_summary`, `summary` (of `Result` and `SweepResult`), `SweepResult` (`results`, `alternatives_table`, `save`), `read_choices`, `write_choices`, `apply_choices`; `RunTables`, `with_choices`, `flagged_points`, `load_saved` | `results.py` (rest), `gui/runs.py` (data part) |
| `src/sweep.jl` | `run_sweep`, `sweep` (`Model.sweep`), `preflight`, `check` (`Model.check`), `preflight_messages`, `point_text` | `sweep.py`, `Model.sweep`, `Model.check` |
| `src/queue.jl` | `solve_chain`, `default_workers`, `Run`, `WorkerPool` (threads: `submit`, `cancel`, `set_workers`, `restart`, `status`, `shutdown`, `on_point`/`on_done`), `sweep_result`, `live_tables`, `load_model_from_args`, `CASE_TEMPLATE` | `runner.py` |
| `src/cli.jl` | The command line: `main(args)` (`@main`, so `julia -m HTGEquilibrium ...`), `parse_cli` (argparse's rules and messages), `cmd_run`, `cmd_check`, `cmd_species`, `cmd_new`; `plot` and `gui` exit with a message | `cli.py` |
| `templates/case_template.toml` | The case template of `htg new` (and the pool's warm-up), a copy of `src/htg/templates/case_template.toml` (identical, tested) | `templates/case_template.toml` |
| `../htg`, `../htg.cmd` | Launchers (Linux/macOS, Windows): `julia --project=<this package> --threads=auto -m HTGEquilibrium ARGS...` | the `htg` console script |
| `test/test_cli.jl`, `test/dump_cli.py`, `test/data/python_cli/` | The command line against Python's output (~70 commands, `to_string` of hand-made tables, an `htg run` and the folder it wrote; 144 kB, committed) | `tests/test_user_species.py` (the CLI flags) |
| `test/test_results.jl`, `test/dump_results.py`, `test/data/python_results/` | Tables and files against a run saved by the Python app (38 kB, committed): the same CSV files and choices.toml, Python's files read as pandas reads them, round trips | `results.py`, `gui/runs.py`, `tests/test_multiple_results.py` (choices) |
| `test/test_sweep.jl` | Grids and chains, preflight messages, warm starts, failed points, solution checks, threads | `tests/test_grid.py`, `sweep.py` |
| `test/test_queue.jl` | The pool with a stand-in chain solver (scheduling, cancel, errors, status, restart, callbacks), real solves | `tests/test_runner.py` |
| `test/compare_cases.jl` | Manual: runs `cases/*.toml` in both versions and compares the saved results | – |

**Python is the reference, not bit for bit** (owner, PR #7: "8e-14 is
considered close"). Arithmetic is plain Julia (`sum`, `^`); floats are compared
with Python within the Python reference tests' tolerances, else rtol 1e-12.
What must match exactly is behaviour: strings and messages, Int vs Float, key
and list order, which values are read and kept. Python's `dict` keeps insertion
order; use `OrderedDict` where order matters (element lists, species lists).

Units follow the Python code: `water_properties` and `hkf_ion_g0` take T in K
and P in Pa and return rho in kg/m3, QBorn and beta in 1/bar, mu0 in J/mol on
the NASA (Hf − T·S) reference. The band functions are for the GUI and take
°C and bar.

```julia
using HTGEquilibrium
w = water_properties(SUPCRT92Water(), 573.15, 240e5)    # w.rho, w.epsilon, ...
mu0 = hkf_ion_g0(SUPCRT92Water(), records, 573.15, 240e5)
supcrt92_band_check([370.0, 378.7, 380.0], 240.0)        # :ok / :in_band / :near_band per T
```

## Build and test

Needs Julia 1.12.7 and gfortran (`apt install gfortran`). From the repository root:

```
julia --project=julia/HTGEquilibrium -e 'using Pkg; Pkg.instantiate()'
julia julia/HTGEquilibrium/deps/build.jl          # or Pkg.build("HTGEquilibrium")
JULIA_NUM_THREADS=4 julia --project=julia/HTGEquilibrium -e 'using Pkg; Pkg.test()'
```

`Manifest.toml` is committed; its versions match `src/htg/juliapkg.json`
(Clapeyron 0.6.30, ForwardDiff 1.5.0, JSON 1.13.4, Roots 3.0.14, ...). The
test environment is `test/Project.toml`, a workspace of this package, so the
same Manifest covers it.

The tests take about 7.5 minutes (5801 tests, 4 threads; the result, sweep,
queue and command-line tests add well under a minute, as the solver is already
compiled); `test_phases.jl` and the
solves are most of it, almost all compilation of Clapeyron (PC-SAFT and
ePC-SAFT, then again for the Duals of the exact Hessian, for each new set of
association sites) and Ipopt. `Pkg.test()` does not pass `-t` on to the test
process; set `JULIA_NUM_THREADS` for the threaded tests.

By default the solver tests solve a representative subset of the reference
points (11 of 21, see `DEFAULT_SOLVES` in `test_equilibrium.jl`); set
`HTG_FULL_TESTS=1` to solve all of them (about 7.5 minutes, 6023 tests, also the CLI run on 2 threads; the extra solves are
fast once compiled):

```
HTG_FULL_TESTS=1 JULIA_NUM_THREADS=4 julia --project=julia/HTGEquilibrium -e 'using Pkg; Pkg.test()'
```

The example cases in both versions (manual, a few minutes; outputs go to a
temporary folder and are not committed; see "Result tables, sweeps and the
worker pool"):

```
JULIA_NUM_THREADS=4 julia --project=julia/HTGEquilibrium julia/HTGEquilibrium/test/compare_cases.jl [--out DIR]
```

The Python fixture of the result tests (`test/data/python_results/`, committed)
is written by `uv run python julia/HTGEquilibrium/test/dump_results.py`.

Optional wider comparison with pychnosz (785 points: 300–600 °C × 200–300 bar,
100–390 °C at 240 bar, and 360–420 °C every 0.5 °C at five pressures). The
output is generated and ignored by git; `runtests.jl` uses it when present:

```
uv run python julia/HTGEquilibrium/test/dump_pychnosz.py --grid
```

### The Fortran library

- **Source.** `vendor/supcrt92/H2O92D.f` is `H2O92D.f.orig` from pychnosz 1.5.4
  (`pychnosz/data/extdata/src/`, as installed in the app's `.venv`). pychnosz's
  README there describes it as the SUPCRT92 file without modification
  (timestamp 15 Nov 1991; J. W. Johnson, LLNL). It is copied byte for byte
  (sha256 `885cf390…cbc3e46`, checked by the tests). It has no licence header;
  the owner accepted this in PR #1.
- **Build.** `gfortran -O2 -fPIC -shared -ffp-contract=off -std=legacy`.
  `-ffp-contract=off` keeps results independent of FMA hardware; `-std=legacy`
  only silences warnings about Fortran 77 constructs. Set `FC` for another
  compiler (the symbol names assume gfortran's `h2o92_`). Set
  `HTG_H2O92_LIB` to load a library from elsewhere. The build output
  (`deps/build/`) is not committed.
- **Call.** `ccall` on `SUBROUTINE H2O92(specs, states, props, error)` with
  the specs pychnosz and CHNOSZ use: `(2, 2, 2, 5, 1, 0, 2, 1, 4, 0)` (°C,
  g/cm3, bar, single phase, T–P input, LVS near the critical point, JN91
  dielectric constant). No wrapper file is needed. The LVS limits used by the
  band check are read from the library's `COMMON /coefs/`.
- **Results.** This library gives the same bits as pychnosz's own
  `h2o92.so` (compared directly at the reference states and in the band
  region). Against the golden data the deviation is zero, see below.
- Outside H2O92's validity range `water_properties` returns all NaN, as
  pychnosz does (the app's `rho >= 350` check then fails, as in Python).
  pychnosz also caches results by T and P rounded to 1e-6; that cache is not
  copied, which matters only for two calls closer than 1e-6 K or bar.
- **Cross-platform builds** (Windows and macOS, ideally a JLL through
  BinaryBuilder/Yggdrasil so users need no compiler) are an **M5 packaging
  task**. Today the library links against the system libgfortran.

### Thread safety

H2O92 is **not reentrant**: all its routines keep state in COMMON blocks
(`/units/`, `/satur/`, `/therm/`, `/param/`, ...) and `SAVE` their locals.
Every call therefore goes through one global lock, `H2O92_LOCK`, so calling
`water_properties(SUPCRT92Water(), ...)`, `hkf_ion_g0` or the band functions
from several threads is safe; the Fortran part runs one call at a time
(~25 µs per call, so this is not a bottleneck next to the solver).

Without the lock it is not just wrong but fatal: a test that called the
library from 4 threads without the lock ended the Julia process within
seconds, through the Fortran `STOP` in `power()` ("neg base \*\* real exp is
complex"), because another thread had overwritten the shared state. With the
lock, the threaded test in `runtests.jl` reproduces the serial results bit for
bit. Results also do not depend on the order of earlier calls (tested).

The same `STOP` (and a second one for `0 ** negative`) would end the process
if H2O92 ever reached it from valid input. It was not reached on any grid
tested (100–600 °C × 150–350 bar, about 1.3 million calls).

`IAPWS95Water` shares one immutable Clapeyron model between threads; the
threaded HKF test gives identical results. Clapeyron's general thread safety
is the subject of M1.

## Golden tests and results

`runtests.jl` (9201 tests, ~17 s including the IAPWS-95 checks after
compilation):

- **thermo.json** (`tests/reference/data/thermo.json`, the three electrolyte
  states 200 °C/240 bar, 300 °C/240 bar, 350 °C/250 bar, both selections):
  the 8 water values with the Python `same` rule (rtol 1e-10, atol 1e-12) and
  `mu0_ion` with `close` (rtol 1e-10, atol 1e-12). **Worst deviation: 0** (bit
  for bit) for water and for ion G.
- **pychnosz grid** (when generated, 785 points): water identical bit for bit;
  ion G worst relative deviation 1.2e-15. HKF follows pychnosz's arithmetic term
  by term, including the parameter scaling round trip of `hkf_ion_g0`.
- Fortran source checksum, NaN outside the validity range, call-order
  independence, threaded calls.

### IAPWS-95

`IAPWS95Water()` is a different equation of state from SUPCRT92's HGK/LVS, so
it is **not** compared with SUPCRT92 or the golden data (the HKF spike found
1e-4 to 6e-3 relative differences in density and up to ~500 J/mol in ion G
near the critical point). Its tests check only that it is self-consistent:
the density equals `Clapeyron.volume` and gives back P; epsilon is JN91 of
that density; alpha, beta, dalpha/dT, Q, Y and X agree with central finite
differences of the model itself (rtol 1e-6 for first, 1e-5 for second
derivatives) at 25 °C/1 bar and 200–450 °C/240–250 bar; HKF with IAPWS-95 is
finite and uses IAPWS-95 water both at (T, P) and at 25 °C, 1 bar.

## SUPCRT92 failure bands (for the GUI warning)

SUPCRT92's LVS density solver (`denLVS`: 20 Newton steps, density clamped to
the LVS limits 199.775 and 420.04 kg/m3) does not converge in bands about
0.5 °C wide on each side of the pseudo-critical line. It then returns a clamp
density without an error flag, with the state up to a few bar off the requested
pressure. At 240 bar:

| band | density returned |
|---|---|
| 378.510–379.016 °C | 420.04 kg/m3 (inside the app's electrolyte range) |
| 386.419–386.948 °C | 199.775 kg/m3 |

**Detection** works at any pressure, without a stored table:
1. call H2O92 at (T, P) and take its density D;
2. call H2O92 at (T, D) to get the pressure of that density on the same
   equation of state;
3. the point has failed if D is exactly a clamp value, or if the implied
   density error beta·|P(T, D) − P| is above 1e-6.

A scan of 365–425 °C every 0.01 °C at 150–350 bar (every 5 bar, ~250,000
points) gave a clean separation. Of the ~43,000 LVS points that did not return
a clamp value, 99 % have an implied density error below 1e-10, all but 21 below
1e-7, 2 lie between 1e-7 and 1e-6, and 19 lie between 1e-6 and 3e-2 (these are
the non-converged ones). Besides the clamp bands there are a few
isolated non-converged points that end without reaching a clamp value
(e.g. 373.74–373.75 °C at 220 bar, density off by ~3 %); the second test
catches those. A scan of 100–600 °C every 0.05 °C at the same pressures found
no failures outside the LVS region (369.85–421.85 °C), where HGK converges to
better than 1e-6 in pressure.

**API** (in °C and bar):
- `supcrt92_failed(T_C, P_bar)`: the point test.
- `supcrt92_failure_bands(P_bar)`: all bands at that pressure. It applies the
  point test every 0.01 °C over the LVS temperature range and refines each
  edge by bisection to 1e-4 °C (~0.3 s per pressure; cached; thread safe).
- `supcrt92_band_check(T_C, P_bar; margin_C=0.5)`: for a sweep. `P_bar` is
  one pressure or one per temperature. Returns per temperature
  `(T_C, P_bar, status, band)` with status `:in_band`, `:near_band` (within
  `margin_C` of a band) or `:ok`. Every temperature is also tested directly,
  so a failure narrower than the scan step is still reported.

The tests check the 240 bar bands against the spike's 0.05 °C scan (and the
spike's bands at 225, 250, 275 and 300 bar), the statuses of a sample sweep,
and that along 100–600 °C at 240 bar the point test fails only inside bands.

## Case files: reading

`load_case(path)` reads a case file in format 2, an old (format 1) file, or a
mixture, with the same rules and messages as Python's `htg.load_case`:

```julia
case = load_case("cases/manure_gas_min.toml")
grid(case), chains(case), temperatures(case)        # sweep points
point = at(case, 300.0)                             # one point (single pressure)
other = replace(case; starts=3)                     # Python: case.replace(starts=3)
dry_element_wt_percent(case.feed, nothing)          # elemental feed: dry wt% used by the model
feed = resolved(case.feed, molar_mass_function(species_elements))   # species-amount feeds
carbon_bound_messages(case, carbon_solids); effective_carbon_bounds(case.carbon_bounds, true)
Units.display(240, "MPa")                           # "24"
```

The structs keep Python's field names and defaults (`CASE_FIELDS`); fields
hold what the file gave (Int or Float64, `OrderedDict`s in file order) until
the constructor normalises the ones the calculation reads, as `__post_init__`
does. Old keys (`gas_carbon_fraction`, `gas_carbon_min_fraction`,
`solid_carbon_fraction`, `phase_carbon_fractions`, `multistart`, `refine`, the
old `[feed]` keys) are read and moved into the new ones on creation.

**Dropped solvers** (owner, PR #10). The Python solvers `"slsqp-fd"` and
`"trust-constr"` and the `[solver]` keys `fd_ftol` and `fd_eps` (settings of
`slsqp-fd`) are gone: `SOLVERS` is `("ipopt", "slsqp")` and `Case` has no
`fd_ftol`/`fd_eps` fields. Old cases still open: a dropped solver is read as
`"ipopt"` (also `Case(...; solver="slsqp-fd")` and `replace`), the two keys are
ignored, and the Julia-only field `case.legacy_notes` says what was changed, for
the app to show, e.g. "The solver 'slsqp-fd' is no longer available; this case
now uses 'ipopt' (the default)." and "[solver] fd_ftol and fd_eps were settings
of the removed solver 'slsqp-fd'; ignored." Writing such a case writes
`ipopt` (the default, so no `method` line) and no fd keys. `legacy_notes` is
not written and not compared by `==`. Files concerned: no case file names
either solver; `cases/yakaboylu2013_fig7_8.toml` (`fd_ftol`, `fd_eps`) and
`cases/cellulose_matlab.toml` (`fd_ftol`), and their old copies in
`tests/data/legacy_cases/`, carry the keys next to `method = "ipopt"`; the
template (`src/htg/templates/case_template.toml`) and its old copy list the
two solvers in a comment only.

**Comparison with Python.** `uv run python julia/HTGEquilibrium/test/dump_cases.py`
writes `test/data/cases_python.json` (~0.6 MB, ignored by git; without it the
comparisons are skipped with an info message). The tests rebuild the same dump
in Julia and compare it field by field: dict key order, Int vs Float, strings
and error messages exactly; floats within 1e-12 relative (or 1e-13 absolute
where Python's value is a cancellation residue such as `100 - sum(wt%)`, 0.0 in
Python and ~1e-14 in Julia). It covers
- the 12 case files: 5 in `cases/`, the 6 old files in
  `tests/data/legacy_cases`, the template. For each: all Case fields,
  temperatures, pressures, grid, chains, `at`, carbon-bound messages, the feed
  as model input (resolved with the registry's element data, which the JSON
  carries for the species the cases name), `feed_fractions`, the feed summary,
  `legacy_to_v2` and the case converted from it;
- 119 variants of the template and old files (one setting or one mistake each,
  66 of them load errors, messages compared);
- 21 [feed] tables and 10 old feeds through the feed functions;
- `format_number`, conversions, `display`, `repr`, `round`, `refine`.

Result: **everything equal; 10519 floats compared, 4 not identical: largest
relative deviation 1.8e-15, one fraction 1.4e-16 where Python has 0.0. The
text `to_toml` writes equals Python's exactly for all 76 cases that can be
written (files, variants, converted old files; 5 more give the same error).**
The Julia-only `water_model` field is left out of the comparison (Python has no
such key; every Python case has the default), and so is `legacy_notes`. The
Python app still has the dropped solvers (above), so its dump (`dump_cases.py`,
unchanged) is mapped on the Julia side before the comparison
(`without_dropped_solvers` in `test/python_cases.jl`): Case fields without
`fd_ftol`/`fd_eps` (a dropped solver would become `"ipopt"`), `to_toml` texts
without their `fd_` lines, and the messages that list the solvers or the
allowed `[solver]` keys without the dropped ones.

**Differences from Python** (deliberate; everything else behaves the same):
- Errors: `CaseError` (a ValueError raised in case.py), `FeedError`,
  `LegacyConversionError` and `Units.UnitError` subtype the abstract
  `ValueError`, so `e isa ValueError` catches what `except ValueError` does.
  (The data layer throws `ArgumentError`s for Python ValueErrors; moving it
  onto `ValueError` is a follow-up. Code that catches Python ValueErrors from
  both should catch `Union{ValueError,ArgumentError}`.)
  Python TypeErrors (`Feed(**data)` with an unknown key, a missing Case field)
  are `ArgumentError`s with Python's text; `replace` with an unknown field
  (Python AttributeError) is an `ArgumentError`.
- Malformed structure that Python reports through an incidental TypeError or
  AttributeError (a section that is not a table, a phase list that is a
  number, a `points` entry that is not a list) gives a readable `CaseError`
  or `ArgumentError` instead; the text differs from Python's.
- TOML syntax errors are `TOML.ParserError`s with the standard library's
  text, not tomllib's.
- `starts = nan` / `inf` and `multistart = nan` / `inf`: a `CaseError` with
  Python's text (Python raises ValueError / OverflowError from `int()`).
- `[aqueous] water_model` (Julia only, see below).
- A molar mass of zero gives `DivideError` (Python ZeroDivisionError).
- `carbon_bound_messages` sums the fixed shares in `PHASES` order (Python
  iterates a set; the sum is compared with an absolute tolerance of 1e-10).
- `registry_elements(registry)` works on the Julia `Registry`;
  `species_elements_function(table)` does the same for any {name => elements}
  table (the tests check `registry_elements(load_model().registry)` against
  the element data Python's registry gives).
- `Units` is a submodule (`Units.display` etc., as `units.display` in Python);
  its conversions broadcast over arrays. `replace(case; ...)` extends
  `Base.replace`.
- `save(case, path)` is not exported (`HTGEquilibrium.save`; the name is
  common). The GUI's `CaseDocument` (tomlkit, keeps comments) is not ported;
  see "Case files: writing".

**Python behaviour kept on purpose.** `round(x, 8)` (applied to every
temperature and pressure; equal values are then merged) → `py_round`, the
correctly rounded decimal; `numpy.arange` of `refine` → `_np_arange` with
NumPy's length rule; `Decimal(repr(x)) + 1` in the legacy feed conversion →
`_plus_one` (gives `1.2`, not `1.2000000000000002`, for the file). Julia's
`TOML.parse` returns plain `Dict`s, so `parse_toml` takes the values from
`TOML.parse` and the key order from a small scanner of the text. Julia's
Printf `%g` chooses between fixed and exponent notation before rounding
(`@sprintf("%.6g", 999999.5)` is `"1000000"`; C and Python give `"1e+06"`);
`py_format_g` follows C for display text, and `formulas.jl` uses it too.

**Water model** (owner, PR #7: a button in the case editor). In `[aqueous]`,
next to the aqueous `model`:

```toml
[aqueous]
model = "electrolyte"
water_model = "SUPCRT92"   # or "IAPWS-95"; absent means "SUPCRT92"
```

`case.water_model` holds the name (default `"SUPCRT92"`, so every existing file
keeps its meaning); `validate` refuses other values ("water_model must be one
of ['SUPCRT92', 'IAPWS-95'], not 'IAPWS95'."). `water_model_of(case)` gives
`SUPCRT92Water()` or `IAPWS95Water()` for `water_properties` and `hkf_ion_g0`
(`WATER_MODELS` maps the names). Python has no such key: it always uses
SUPCRT92, and its `load_case` would refuse the key as unknown.

## Case files: writing

**Decision (owner, PR #8): write a fresh file; comments are not kept.** Users
will mostly edit cases in the GUI, not the files. Keeping comments or other
context (e.g. in a separate file) may come later.

`to_toml(case)` / `dump_case(case)` give the text, `HTGEquilibrium.save(case,
path)` writes it (Python: `case.to_toml()`, `dump_case`, `case.save`). As in
Python: format 2, `format` and `name` first, then `description`,
`atomic_weights`, `[conditions]`, `[feed]` (+ `[feed.species]`,
`[feed.salts.X]`), `[phases]`, `[aqueous]`, `[constraints]`
(+ `[constraints.carbon]`), `[initial_guess]`, `[solver]`, `[sweep]`,
`[comparison]`; only values that differ from the defaults (plus the required
ones); numbers as Python's `repr` (`1.0`, `1e-06`), Int stays Int; old (format
1) settings written in the new form (`refine` merged into `temperature_C`, an
old feed converted to format 2 where possible); `comparison.file` as an
absolute path. `water_model` (Julia only) appears in `[aqueous]` only when it
is not "SUPCRT92", so other files are what Python writes.

Tests: every file in `cases/` and the template round-trips (load, save, load
gives an equal `Case`); old files give format 2 with the same model numbers;
the written text has no comments, while the template it came from has many;
with `test/data/cases_python.json`, the text equals Python's `to_toml` for
every file, variant and converted old file.

## Data layer (databases, registry, standard states)

`load_model(dir)` is Python's `Model.load(dir)`: it reads the four databases,
`species.toml` (with an optional `user_data` folder layered on top), discovers
the species and returns a `Model` (`databases`, `registry`, `user_info`,
`load_args`). It takes ~0.1 s after compilation. Nothing is cached globally; a
loaded `Model` is only read afterwards, so one model can be shared by all
threads and all cases (tested: threaded and serial reads give identical results).

**Water model.** The `Model` holds no water model. The water model
(`SUPCRT92Water()` or `IAPWS95Water()`) is a per-case choice: the case key
`[aqueous] water_model`, `water_model_of(case)` (see "Case files"), chosen
with a button in the case editor (owner, PR #7; the GUI is a later step). Every calculation that needs water properties takes
it as its first argument (`water_properties(water_model, T, P)`,
`hkf_ion_g0(water_model, model, ions, T, P)`), so a case uses one model
throughout and never mixes them (owner decision 3).

```julia
model = load_model()                                  # data/ (or $HTG_DATA_DIR)
model = load_model("tests/reference/databases"; nasa_names=Dict(), user_data=resolve_user_dir())
species(model; include_unavailable=true)              # Vector{CatalogueRow}
resolve_parameters(model.registry, model.databases, 623.15, ["water", "methane"])
nasa_g0(model, ["water"], 623.15); solid_g0(model, ["Calcite"], 623.15, 240e5)
hkf_ion_g0(SUPCRT92Water(), model, ["Na+", "HCO3-"], 573.15, 240e5)   # the case's water model
```

**Results.** `catalogue.json` (both reference models: the catalogue of 1901 /
1900 rows, the four registries with their order, the ion ePC-SAFT and HKF
tables, the mineral inventory, the product candidates) is reproduced fully
with Python's `same` rule, except the intended gas candidates (below). In
`thermo.json` the molecular g° (235 species), the solid g° (169, including the
Python error texts for minerals outside their range), the SAFT parameters
(26), and mu0 of the fluid phases and of the electrolyte phase (ions and
neutrals) agree at all six states within the Python tests' tolerances; worst
relative deviations: molecular g° 8.4e-14, fluid mu0 3.1e-16, solid g°
4.1e-16, SAFT parameters, mu0 of ions and neutrals 0. A user species folder
written by the Python app (`test/dump_user_species.py`) loads to the same
catalogue rows, registries, user_info, parameters and g° (rtol 1e-12).

Floats are plain Julia arithmetic (`^`, `sum`), not Python's bits (owner,
PR #7): tests compare them with the tolerances of the Python reference tests
(`tests/helpers.py`: rtol 1e-10, atol 1e-12) or rtol 1e-12 elsewhere; strings,
integers and orders are compared exactly.

**Gas product candidates (a change from Python, owner PR #7).**
`product_candidates` lists a molecule as a gas candidate when its NASA-7
record is a gas record: phase letter (header column 45) "G" **or blank**. In
`Thermdat_WasteMech.txt` the gas-phase mechanism species (CH4, H2, N2O, ...;
187 of 229 records) leave that column blank, only TARN has "G", and the 41
solids have "S" (CHEMKIN treats such a gas-mechanism file the same way). Python
required "G", so it never listed a gas candidate. With the default salt
elements at 400 °C this adds 7 gas rows before the 43 solids: chlorine (ready)
and NOCL, CLCO, CL, CLO, HCL, HOCL (NASA records without SAFT parameters; not
ready, with the reason).

**Differences from Python** (names are kept; the rest follows Julia habits):
- Methods become functions: `Databases.load` → `load_databases`,
  `Model.load` → `load_model`, `UserDatabase.load` → `load_user_database`,
  `registry.all_molecular` → `all_molecular(registry)`, `model.species()` →
  `species(model)`, and functions that change their argument end in `!`
  (`discover!`, `apply_user_layer!`, `hide_species!`).
- Exceptions: Python `ValueError`/`TypeError`/`FileNotFoundError` are
  `ArgumentError`s, `KeyError`s are `KeyError(message)`; the messages are
  Python's (including texts such as "could not convert string to float: ''"
  that end up in the catalogue).
- pandas tables: the SAFT CSVs are `DataTable`s that reproduce pandas'
  reading (default NA strings, numeric columns as Float64 with NaN, text
  columns as String/`missing`, the index labels, user rows first with negative
  labels). The catalogue is a `Vector{CatalogueRow}`; the ion tables, parameter
  sets, registry entries and inventory rows are `OrderedDict{String,Any}` with
  Python's keys.
- `supcrt92.xml` is parsed with EzXML into an `XMLNode` tree with
  ElementTree's semantics (`findtext`, `.//tag`, `iter`).
- TOML: Julia's TOML stdlib returns unordered `Dict`s, but the registry order
  comes from `species.toml`; `read_toml` (toml_ordered.jl, shared with the
  case files) rebuilds the tables as `OrderedDict`s in file order (values from
  the stdlib parser; each table of an array of tables in its own order, as tomllib).
- Gas product candidates (above); their Reason is the plain error text.
- Not ported (later steps): `species_catalogue.py` (Species tab table, M4),
  editing/importing user species (M4). (`Model.solve`, `Model.sweep` and
  `Model.check` are `solve`, `sweep` and `check`.)

**For the phases step.** `resolve_parameters` gives the parameter sets
(`m`, `sigma`, `epsilon`, `epsilonAB`, `kappaAB`, `association_scheme`, `Mw`)
that `build_pcsaft_model` needs; components are `all_molecular(registry)[name]["formula"]`
with `registry.binary_interactions`. The electrolyte phase needs
`aqueous_reference_neutrals/ions`, `probe_weights`, `ion_charge`,
`registry.ion` (formula, elements), `ion_epcsaft_table` (m, sigma, epsilon,
dielectric_constant), `hkf_ion_g0(water_model, model, ions, T, P)`,
`nasa_g0(model, neutrals, T)` and `water_properties(water_model, T, P)` for the
validity check and the dielectric constant (Python calls `supcrt_water` there;
Julia takes the case's water model). The phase parts of `thermo.json` are
tested in `test_phases.jl` (see "Phase models").

## Phase models (`phases.jl`)

Port of `phases.py` with Python's names, formulas and messages. Clapeyron is
used directly, with the arguments the Python app passes through juliacall
(`src/htg/backends.py`): `PCSAFT(components; userlocations)` for fluid
phases; for the electrolyte `ePCSAFT(components, charges, short_range.idealmodel,
short_range, hsdDH(components, ConstRSP(eps_r), String[]), String[])` with
`short_range = PCSAFT(components; userlocations, verbose=false)`. While solving,
the only Clapeyron call is `fugacity_coefficient(model, P, T, n; phase)`, as in
Python. The M1 spike (`julia/spikes/m1_threads/replay.jl`) was the starting
point; the phases now take their data from the registry and standard states.

```julia
model = load_model()
T, P = 573.15, 240e5
names = ["water", "hydrogen", "carbon dioxide", "methane"]
resolved = resolve_parameters(model.registry, model.databases, T, names)
eos = build_pcsaft_model(names, resolved, model.registry.binary_interactions)
gas = FluidPhase("gas", names, nasa_g0(model, names, T), eos, T, P, "vapour")
G, mu = gibbs_and_mu(gas, n); lnphi(gas, n); hessian(gas, n)
solid = SolidPhase(minerals, solid_g0(model, minerals, T, P))
# resolved must also hold registry.aqueous_reference_neutrals
aq = ElectrolytePhase(water_model_of(case), model, selected, T, P, resolved, ATOMIC_WEIGHTS)
gibbs(aq, n); gibbs_and_mu(aq, n); aq.lnphi_infinity; HTGEquilibrium.phi(aq, full)
```

**API** (`Phase` is the abstract type of the three):

| | `SolidPhase` | `FluidPhase` | `ElectrolytePhase` |
|---|---|---|---|
| `gibbs(p, n)`, `gibbs_and_mu(p, n)` | ✓ | ✓ (empty phase: G = 0 and the pure-component mu, cached) | ✓ (no water: 0, zeros) |
| `mu_at(p, n)`, `lnphi(p, n)`, `hessian(p, n)` | – | ✓ | ✓ |
| `phi(p, full)`, `gibbs_full(p, full)`, `expand(p, n)` (not exported) | – | – | ✓ |

Fields: `species` (all); `name`, `mu0`, `eos_model`, `T`, `P`, `root` (fluid);
`water_model`, `T`, `P`, `species`, `neutrals`, `ions`, `all_species`,
`positions` (1-based indices of `species` in `all_species`), `z`,
`probe_weights`, `water` (a `WaterProperties`), `components`, `charges`,
`short_range`, `eos_model`, `mu0_ion`, `mu0_neutral`, `lnphi_infinity`
(electrolyte). The functions accept any real element type, so ForwardDiff can
differentiate through them.

**Errors.** `EOSFailure <: ValueError` (Python: `class EOSFailure(ValueError)`):
no valid fugacity coefficients; on the optimizer path (`gibbs_and_mu`,
`lnphi`, `hessian` of the electrolyte) every probe failure becomes an
`EOSFailure`, as in Python. `PhaseError <: ValueError` stands for the plain
Python `ValueError`s of `phases.py` (same messages); missing parameter sets are
`KeyError(message)` as in the data layer.

**Water model** (owner decision, PR #7). `ElectrolytePhase` takes the case's
water model as its first argument (`water_model_of(case)`; Python has no
argument and calls `supcrt_water`). That one model gives the validity check
(`rho >= 350 kg/m3`, `epsilon > 1`), the dielectric constant of the
Debye-Hückel term and, through `hkf_ion_g0(water_model, model, ...)`, the ion
g°; `phase.water_model` and `phase.water` record it. Tested: with
`IAPWS95Water()` the phase's water properties are IAPWS-95's, mu0 of the ions
equals `hkf_ion_g0(IAPWS95Water(), ...)`, and G differs from SUPCRT92's.

**Threads** (M1 decision). No global state. Clapeyron models are only read,
so one model may be shared by several phase objects and threads. An
`ElectrolytePhase` is not changed after construction. The only cache,
`FluidPhase`'s pure-component mu of an empty phase (per phase object, as in
Python), is filled under the phase's own lock. The intended pattern: each task
builds its own phase objects; it may reuse a fluid phase's Clapeyron model.
Building an `ElectrolytePhase` (its two Clapeyron models and the
infinite-dilution reference) takes ~0.5 ms after compilation. Tested: 24 tasks
on 4 threads (the gas phase of all six `thermo.json` states with one shared
Clapeyron model per state, the electrolyte phases of three states built per
task; G, mu, ln phi, Hessians, empty phases) give results identical bit for
bit to the serial run; 16 tasks filling one shared phase's cache get the same
values.

**Results against `thermo.json`** (tolerances of `tests/test_reference_thermo.py`:
G rel 1e-10 / abs 1e-8, ln phi and `lnphi_infinity` rtol 1e-9 / atol 1e-11;
6 states × 3 fluid phases × 5 points, 3 states × 2 electrolyte selections × 5
points). All pass; worst relative deviations: fluid G 3.1e-16, fluid ln phi
2.2e-16, electrolyte G 6.0e-16, electrolyte ln phi 3.1e-13 (on a small value;
absolute < 1e-14), `lnphi_infinity` 2.4e-14. The 7 reference EOS failures
(fluid phases at 200, 350 and 375 °C) raise `EOSFailure` with Python's message
(the composition list within 1e-12). Ported from `tests/test_phases.py`: mu =
dG/dn by central differences (and for the Na+ + HCO3- pair), Euler
`G = n·mu`, homogeneity of degree one.

**Derivatives** (owner, PR #9: exact derivatives). mu comes from the
fugacity coefficients, as in Python; it is the gradient of G (central
differences, and ForwardDiff to 4e-15 in the first version of this port).
`hessian(phase, n)` is exact: the ForwardDiff Jacobian of `mu_at` through
Clapeyron, symmetrised. Python used an analytic ideal part plus forward
differences of ln phi; that version (`_hessian`) is gone. Zeros for an empty
phase (fluid: sum(n) <= 1e-13; electrolyte: water <= 1e-9 mol).

For the electrolyte, a single ion cannot be varied alone: the EOS probe makes
every composition electroneutral (smoothly, through `hypot(imbalance, 2e-14)`),
so a single-ion column of the Jacobian is the derivative along that ion with
half of its charge moved to the first counter-ions. Along electroneutral
directions Z these terms cancel, so J·Z is exact; the symmetrised H keeps
Z'HZ, which is what the optimizer needs (its balances keep the aqueous phase
electroneutral). Python's forward differences got Z'HZ wrong by ~16 % along
the Na+ + HCO3- pair (each one-ion step was neutralised by the first counter-ion,
so the pair was counted twice).

`_exact_hessian` seeds ForwardDiff Duals itself, 8 columns per evaluation
(`HESSIAN_CHUNK`), with one tag per phase type, so every phase size uses the
same Dual type. Compilation is the cost: the first Hessian of a session takes
~25 s for PC-SAFT and ~40 s for ePC-SAFT, and a Clapeyron code path met for the
first time with Duals (e.g. more association sites) adds ~20 s once. The app
should precompile or warm this up at start (M5). Tests (`test_phases.jl`, 300 °C, 240 bar):

| check | worst relative deviation |
|---|---|
| fluid Hessian vs central differences of mu (3 phases × 2 points, column by column) | 5.7e-9 |
| fluid H·n = 0 (Gibbs-Duhem), relative to max(|H||n|) | 1.3e-13 |
| electrolyte Z'HZ vs central differences of Z'mu (5 and 16 selected species; neutrals, Na+ + HCO3-, 2 Na+ + CO3^2-, K+ + Cl-, Mg2+ + 2 Cl-, Cs+ + Li+ + CO3^2-) | 1.2e-8 |
| electrolyte (H Z) rows of the neutral species vs central differences | 3.7e-9 |

Symmetry is exact (H == H'); the electrolyte also passes Z'H n = 0.

Cost per call after compilation (one point, 1 thread):

| phase | `gibbs_and_mu` | `hessian` |
|---|---|---|
| gas, 5 species | 7 µs | 17 µs |
| gas, 9 species | 18 µs | 123 µs |
| electrolyte, 5 selected (11 components) | 33 µs | 93 µs |
| electrolyte, 16 selected | 49 µs | 286 µs |

For comparison, Python's forward-difference Hessian took 47 µs (gas, 5
species) and 215 µs (electrolyte, 5 selected) in this port; it needs k + 1
evaluations of ln phi, the exact one ceil(k / 8) evaluations with Duals.

**Differences from Python** (deliberate):
- `ElectrolytePhase(water_model, model, selected, T, P, resolved, atomic_weights)`:
  the water model comes first (above). The validity message keeps Python's
  text, "State outside SUPCRT aqueous validity checks", also for IAPWS-95.
- Methods are functions (`gibbs(phase, n)`); `positions` are 1-based;
  `phase.water` is a `WaterProperties` (`water_dict` gives Python's dict);
  Python's `hasattr(model, 'all_species')` is `phase isa ElectrolytePhase`.
- `PhaseError` for Python's plain ValueErrors (above).
- `np.polyfit(x, y, 1)` of the infinite-dilution extrapolation is the
  centred least-squares line (the same line; `lnphi_infinity` within 2.4e-14).
- `SolidPhase`'s `gibbs_and_mu` returns a copy of mu0 (Python returns the
  array itself).
- The pure-component cache of `FluidPhase` is filled under a lock.
- `hessian` is exact (ForwardDiff), not forward differences (owner, PR #9).
- The electrolyte's emptiness checks use the water amount found by name
  (`water_index`, the position of "water" in `species`; 0 when not selected),
  never the first selected species (owner, PR #9; Python's `hessian` tested
  `n[0]`).

**Kept from Python:** an empty gas phase containing water raises `EOSFailure`
where pure water has no vapour root (e.g. 300 °C, 240 bar), because the
pure-component mu is evaluated at pure water.

**For the optimizer and equilibrium steps.** `GibbsProblem` needs, per phase:
`gibbs_and_mu(phase, n)` (the objective and gradient; an `EOSFailure` makes
`fun_and_grad` return NaN), `hessian(phase, n)` (`hessian_blocks`, exact; skip
`SolidPhase`; an `EOSFailure` there means a NaN Hessian, as in Python), `mu_at(phase, n)` (tangent-plane test), `gibbs(phase, n)`
(`phase_G`, the G check) and `phase.species`, `phase.T`. Equilibrium builds the
phases as `equilibrium.py` does: `resolve_parameters` over all non-ion species
(plus the reference neutrals when the aqueous model is electrolyte), then
`ElectrolytePhase(water_model_of(case), model, ...)`, `SolidPhase(names,
solid_g0(model, names, T, P))`, `FluidPhase(phase, names, nasa_g0(model,
names, T), build_pcsaft_model(names, resolved, registry.binary_interactions),
T, P, root)`. Catch `EOSFailure` where Python catches it and `ValueError`
(the abstract type) where Python catches `ValueError`; note that the data
layer still throws `ArgumentError` for Python ValueErrors.

## Solver (`results.jl`, `optimizer.jl`, `equilibrium.jl`)

```julia
model = load_model()                   # the app's data/ folder
r = solve(model, case)                 # Model.solve: one point, or a list point by (T, P)
r = solve(model, case, 450.0, 250.0)   # one point of a sweep case
r.gibbs_J_per_mol_C, r.amounts, r.optimizer["starts"], r.water_model
s = setup(model, case)                 # the Setup (problem, feed, bounds) without solving
```

Everything Python's `equilibrium.solve` does is there: the phase rules
(cut-off, disabled aqueous phase), the feed (elements or species, water,
salts), carbon bounds, the starting points (LP, phase, random, user amounts),
IPOPT with SLSQP as the fallback, multistart and grouping into alternatives,
the polish step and its safeguard, the tangent-plane test, the G check against
the phase models, and the KKT diagnostics. Failures are `SolveError`s whose
message starts with the stage, as in Python (`"feed: ValueError: ..."`;
Python's exception names are kept in the message).

**Dependencies** (exact compat):

| Package | Version | Why |
|---|---|---|
| Ipopt.jl / Ipopt_jll | 1.16.0 / 300.1400.1902 | IPOPT (C API), the default solver |
| MUMPS_seq_jll, SPRAL_jll, OpenBLAS32_jll | 500.900.100, 2025.9.18, 0.3.34 | IPOPT's linear solvers, as in `src/htg/juliapkg.json` |
| NLopt.jl | 1.2.1 | SLSQP (`LD_SLSQP`), also for the tangent-plane test and the polish |
| HiGHS.jl / HiGHS_jll | 1.20.1 / 1.12.0 | The LP starting points; 1.12.0 is the HiGHS in SciPy 1.18.1, which gives the same LP vertices as Python (HiGHS 1.16 picks other vertices of degenerate LPs) |

**Against `solves.json`** (tolerances of `test_reference_solves.py`):

- Setup of all 21 points (initial amounts, feed atoms, C per g, phase G): the
  deviation is 0 except phase G (at most 1.5e-15 relative).
- Points solved with IPOPT in the reference: G within 1e-15, compositions
  within 1e-12 (worst 1.8e-10 of the tolerance). SLSQP (`slsqp_exact`): G
  within 1.6e-13, compositions within 1.3 % of the tolerance.
- Points solved with the original finite-difference SLSQP (`slsqp-fd`, dropped)
  are solved with IPOPT. The stored reference of these points came from a bad
  optimizer and is **a rough guide only** (owner, PR #10): the G must not be
  higher (worst +9e-11 relative). Four points have a lower minimum than the reference
  (`LOWER_THAN_FD`): cellulose_400C (-2.2e-6 relative), full_text_250C
  (-2.6e-9), full_text_gas_min_300C (-6.9e-9), configurable_300C (-1.1e-9).
  For these the compositions are compared with the Python app's own IPOPT
  solution (`test/data/solves_python_ipopt.json`, written by
  `test/dump_solves.py`), which they match (G within 1e-15, amounts within
  1e-10). IMPROVED (full_text_250C_ipopt) and FD_STOPS_EARLY (fig78_380C)
  are kept as in Python.
- Expected failures give the same stage and message; the polish safeguard
  still triggers (stationarity 12 RT, polish used, G within 3e-9 J).
- IPOPT iterations: mostly identical to Python; with the exact Hessian a few
  differ by one or two (fig78_450C_all_solids 23 vs 21,
  cellulose_species_feed_450C 24 vs 25, configurable_300C_ipopt 34 vs 33).

**Threads:** `solve` may be called from several threads (tested: 24 solves, 36 with `HTG_FULL_TESTS`, on
4 threads equal the serial results). The HiGHS calls are under a lock; nothing
else is shared. The tests set `BLAS.set_num_threads(1)`; the app should too
(threads per solve, not inside the linear algebra).

**Differences from Python:**

- **`trust-constr` and `slsqp-fd` are dropped** (owner, PR #10). A case that
  names either is read as `ipopt` with a note (see "Case files: reading");
  `optimize` knows only "ipopt" and "slsqp" ("Unknown optimizer 'slsqp-fd';
  choose from ['ipopt', 'slsqp']").
- SLSQP is NLopt's `LD_SLSQP`, not SciPy's: the same algorithm, other
  iteration counts and messages ("Optimization terminated successfully" /
  "Iteration limit reached" are mapped; others are NLopt's).
- **`max_iterations` for SLSQP** (owner, PR #10): SciPy counted iterations,
  NLopt counts function evaluations (each iteration takes one or more). SLSQP
  stops after `SLSQP_EVALUATIONS_PER_ITERATION * max_iterations` = 10 ×
  `max_iterations` evaluations: 50000 with the default 5000, and a value set in
  the case scales the same way (so it still limits roughly the iterations, as
  in Python). The tangent-plane test (SciPy: 200 iterations) gets 2000
  evaluations. For IPOPT `max_iterations` is IPOPT's `max_iter`, as before.
  `Result.iterations` and the `nit` of a start count NLopt's evaluations.
- Random starts use a seeded `Xoshiro`, not NumPy's stream; the starts differ,
  the seed still makes a run repeatable.
- The IPOPT Hessian is exact (`hessian_blocks`, owner PR #9); Python passed a
  finite-difference Hessian to IPOPT. The IPOPT options are Python's
  (`_JULIA_IPOPT` plus the case's `ipopt_options`).
- Indices are 1-based (`Setup` positions, `group_solutions` members).
- `Result` has a `water_model` field (the case's water model).
- An error other than `EOSFailure` inside an IPOPT callback stops IPOPT and is
  rethrown after the solve (Ipopt.jl does not catch exceptions in callbacks).
- `EquilibriumError <: ValueError` for Python's plain ValueErrors,
  `RuntimeError` for Python's RuntimeErrors.

## Result tables, sweeps and the worker pool (`tables.jl`, `results.jl`, `sweep.jl`, `queue.jl`)

```julia
model = load_model()
case = load_case("cases/manure_gas_min.toml")
check(model, case)                                   # Model.check: preflight table (no EOS)
sw = sweep(model, case)                              # Model.sweep; progress lines as Python prints them
sw = sweep(model, case; workers=4)                   # chains on 4 threads (julia -t 4)
summary(sw); species_table(sw); phase_summary(sw); alternatives_table(sw)   # DataTables
HTGEquilibrium.save(sw, "out/manure_gas_min")        # summary.csv, species.csv, phases.csv, alternatives.csv, case.toml

pool = WorkerPool(4; model)                          # or model_args=Dict("data_dir" => ...)
run = submit(pool, case, "out/manure_gas_min")       # saved there when the run finishes
wait(run); run.status, run.done_points, run.failed_points
status(pool); cancel(pool, run); set_workers(pool, 2); restart(pool, load_model())
shutdown(pool)

tables = load_saved("out/manure_gas_min")            # a saved run of either version, choices applied
write_choices("out/manure_gas_min", Dict((400.0, 240.0) => "B"))
```

**Tables and CSV files.** The tables are the data layer's `DataTable`, with
one Julia type per pandas dtype (Float64 with NaN for float64, Int for int64,
Bool for bool, String/`missing` for text). `DataTable(rows)` infers the types
from a list of dicts as `pd.DataFrame(rows)` does: columns in order of first
appearance, a missing key is NaN, an int column with a missing value becomes
float (so `iterations` is written `29.0` when a sweep has a failed point, as
pandas writes it), a bool column with a missing value becomes object. No new
package was needed: `csv_text`/`to_csv` is a small writer that reproduces
`DataFrame.to_csv(index=False)` (pandas 3.1.1): floats as Python's repr
(`300.0`, `1e-06`, `1.2345678901234568e+17`), NaN and missing as empty
fields, bools as `True`/`False`, text quoted only when it holds a comma, a
quote or a line break (`csv.QUOTE_MINIMAL`; a row that is one empty field is
`""`), `\n` line ends. `read_csv` reads with CSV.jl (as text) and applies
pandas' type inference (default NA strings; int, float, bool or text).

**Compatibility of saved runs (M3).** Same file names, columns, column order,
row order and number formatting as Python:

- test/data/python_results (Python's `SweepResult.save` and `write_choices`
  of hand-made Results with ions, an empty phase, a failed point whose message
  has a comma, quotes and a line break, a removed aqueous phase with a dropped
  bound, a missing `descent_RT`, points with 2 and 3 results): Julia builds the
  same Results and writes the **same bytes** for summary.csv, species.csv,
  phases.csv, alternatives.csv, case.toml and choices.toml.
- Python's files read in Julia (`load_saved`, choices applied, and the tables
  as saved) give Python's columns, dtypes and values (python_read.json); Julia's
  files read back to the tables Julia wrote, and writing them again gives the
  same bytes.
- Every CSV of the 5 example runs written by Python (below), read by
  `read_csv` and written by `csv_text`, gives the original file byte for byte.

**Example cases in both versions** (`test/compare_cases.jl`, manual; see its
header for the rules). Python: `uv run htg run CASE --out DIR --no-plots
--no-user-species`; Julia: all five cases submitted to one `WorkerPool` on 4
threads, which saves each run when it finishes. Result (2026-11-24):

| case | points | match | notes |
|---|---|---|---|
| cellulose_matlab | 8 | 8 | worst G 1.9e-7 of the tolerance (rel 1e-9); amounts 9e-10 of it |
| manure_configurable (continuation) | 31 | 31 | worst G 1.9e-6 of the tolerance |
| manure_gas_min | 16 | 16 | worst G 1.6e-6 of the tolerance |
| manure_ions_minerals | 16 | 16 | worst G 0.34 of the tolerance (3.4e-10 relative); amounts 4.4 % of it |
| yakaboylu2013_fig7_8 | 16 | 16 | 389.6 °C fails in both with the same message (SUPCRT aqueous validity) |

All 87 points agree; no point has a lower G in either version; the
failure message is the same text. File names, columns and rows are equal;
case.toml is equal except that Python still writes `fd_ftol`/`fd_eps` (the
dropped solver's keys) for cellulose_matlab and yakaboylu2013_fig7_8.
IPOPT iteration counts differ at 23 of the 87 points (exact Hessian; usually
by 1–2, at most 20 at manure_ions_minerals 340 °C: 54 vs 34), which is
reported, not compared. Time: Python 60–84 s per case (387 s in total; each
run starts its own Julia, about 40 s of that); Julia 170 s for all five
together on 4 threads, including about 2 minutes of compilation.

**Sweeps.** `run_sweep(model, case; progress=println, workers=1)` and
`sweep` (`Model.sweep`) solve the chains with `solve_chain`, which starts each
point of a chain from the previous solution when the phase selection is the
same (`"solution at 300 °C, 240 bar"`), else from a fresh LP (`"feasible LP
(phase selection changed)"`); failed points are recorded as their message.
With `workers > 1` and more than one chain the chains are solved on a
`WorkerPool` sharing the model. `preflight`/`check` (`Model.check`) and
`preflight_messages` check every point without the equation of state, with
Python's messages.

Tests (`test_sweep.jl`, `test_queue.jl`): the warm-start chain through the
aqueous cut-off (progress lines as Python's; warm starts reach the same minimum
as cold starts, G within 1e-9); every sweep solution passes checks that need no
reference (element balances, amounts >= 0, electroneutral aqueous phase, carbon
shares and phase G add up); **threads**: a sweep on 4 workers equals the serial
sweep exactly (amounts, G, iterations) for 4 continuation chains of 2 points
and for 8 independent points (3 temperatures with `HTG_FULL_TESTS`); the pool's
scheduling (FIFO, one chain per worker, cancel, internal errors, status,
`set_workers`, `restart`, callbacks, warm-up, shutdown) with a stand-in chain
solver, as `tests/test_runner.py`; real solves on 2 workers equal `solve`.

**Worker pool.** One Julia process; each worker is a task on a thread
(`Threads.@spawn :default`), taking chains from one FIFO queue; `workers`
limits the number of chains solved at once (start Julia with enough threads).
User-facing behaviour as Python's `WorkerPool`: `Run` (`id`, `case`,
`output_dir`, `total`, `status` queued/running/cancelling/done/cancelled,
`outcomes`, `log`, times, `run.name`, `run.done_points`, `run.failed_points`,
`run.finished_running`, `wait(run, timeout)`, `sweep_result(run)`),
`submit(pool, case, output_dir)`, `cancel` (queued pieces dropped, running
chains finish), `set_workers`, `status` (per worker: id, state, message, case,
point being solved; queued pieces; target), `shutdown`, `on_point(run, T, P,
outcome, start)` and `on_done(run)`, `default_workers()`. When a run finishes
it is saved to its `output_dir` (if it has one and any point has an outcome),
then `on_done` is called, as in Python. `BLAS.set_num_threads(1)` is set when
the pool starts and in `run_sweep`.

**Differences from Python:**

- Threads instead of processes (M1). The model is loaded once and shared; no
  per-worker Julia start-up. `warm_up` solves the template once in one worker
  (compiled code is shared by all threads).
- `restart(pool, model_args)` (or `restart(pool, model)`) loads the new model
  and uses it for runs submitted from then on, at once; runs already submitted
  finish with their own model (`run.model`). Python made later runs wait for
  new worker processes. Both never mix data within a run. Returns true when no
  run was unfinished.
- A worker cannot crash the way a process can: an error escaping
  `solve_chain` (other than a `SolveError`) fails the rest of that chain with
  "internal error: ..." (Python: the same for an exception in the worker; a
  dead process gave "worker process stopped unexpectedly"). A crash of Julia
  itself (e.g. Fortran `STOP`) would end the whole app (issue #5).
- `on_point`/`on_done` are called from a worker thread (under the pool's lock,
  one at a time); Python called them from its dispatcher thread.
- `shutdown` lets each worker finish the point it is solving (Python waited 5 s,
  then terminated the processes); the rest of its chain is not solved.
- `default_workers()` is min(Julia threads, CPU threads − 1, free memory /
  0.25 GB) (Python: processes of 1.5 GB each).
- The progress line of a parallel sweep says "Solving N points with W workers
  (threads; T available)" (Python: "... in W worker processes (each starts
  Julia; ...)").
- Python methods are functions (`summary(sweep)`, `species_table(r)`,
  `submit(pool, case)`, `HTGEquilibrium.save(sweep, dir)`); `summary` extends
  `Base.summary`.
- Not yet written: the user species run record (`user_species_used.toml` and
  the note at the end of case.toml, `user_species.run_record`), which needs
  the user species editing port (M4). Runs made with a user species folder
  therefore lack it in Julia.
- The water-density check of `preflight` uses the case's water model (Python:
  SUPCRT92 always), as the electrolyte phase does.

## Command line (`cli.jl`)

The commands, arguments, defaults, printed text and exit codes of Python's
`htg` (`src/htg/cli.py`): `run`, `check`, `species`, `new`, with `--data DIR`
(before the command), `--user-data DIR` and `--no-user-species`. From the
repository root:

```
julia/htg species --search methane              # Linux, macOS
julia/htg check cases/manure_gas_min.toml
julia/htg run cases/manure_gas_min.toml --workers 4 --no-user-species
julia/htg new my_case.toml

julia\htg.cmd run cases\manure_gas_min.toml --workers 4     # Windows (cmd or PowerShell)
```

The launchers `julia/htg` (sh) and `julia/htg.cmd` run
`julia --project=<repo>/julia/HTGEquilibrium --threads=auto --startup-file=no -m HTGEquilibrium ARGS...`
(`JULIA_NUM_THREADS` chooses another thread count, `JULIA` another julia
executable). Without a launcher:

```
julia --project=julia/HTGEquilibrium -t auto -m HTGEquilibrium run cases/manure_gas_min.toml --workers 4
julia --project=julia/HTGEquilibrium -e 'using HTGEquilibrium; HTGEquilibrium.main(["check", "cases/manure_gas_min.toml"])'
```

`main(args; out=stdout, err=stderr)` returns the exit code; `@main` makes it the
entry point of `julia -m HTGEquilibrium`, whose exit code it becomes. Each
command starts Julia and loads the package (about 10 s on this machine before
anything is printed); the solver then compiles on the first point (about a
minute), as in Python, where each command started Julia through juliacall.
The package is precompiled on first use after a change (`Precompiling
packages...`, ~5 s).

**Threads (`--workers N`).** `--workers N` now means N worker threads of the
thread-based `WorkerPool` (Python: N worker processes, each starting Julia).
Julia's thread count is fixed at start-up, so the CLI uses min(N, threads Julia
has) and says so when it has fewer ("Julia has 4 thread(s); using 4 worker(s)
instead of 8 (start Julia with more threads: ...)"). The launchers start Julia
with `--threads=auto` (one thread per CPU thread), so N is only capped where
more workers than CPU threads would not help anyway. This was chosen over the
alternatives: the launcher reading `--workers` itself (argparse abbreviations
such as `--work=4` make that fragile) and the CLI restarting Julia with `-t N`
(a second start-up of ~10 s and more code). `--workers 1` (the default) solves
in the main thread, as `model.sweep` does in Python; idle threads cost nothing.

**Pkg app** (Julia 1.12, `[apps.htg]` in Project.toml with
`julia_flags = ["--threads=auto"]`). Optional:
`julia -e 'using Pkg; Pkg.Apps.develop(path="julia/HTGEquilibrium")'` writes a
shim `~/.julia/bin/htg` (Windows: `htg.bat`) that runs this checkout with the
committed Manifest; with `~/.julia/bin` on PATH, `htg run ...` works from any
folder. Tried in a separate depot (shim, arguments and exit codes work). Apps
are new and marked experimental in Pkg 1.12, the shim passes arguments before a
`--` to Julia, and they install into the user's depot, so the launchers stay
the documented way for now; apps are worth adopting in M5 (packaging), when the
package can be installed from a URL or registry with `pkg> app add` and the
Fortran library comes as a JLL.

**Printed tables.** `to_string(table)` (tables.jl) reproduces pandas 3.1.1's
`DataFrame.to_string(index=False)` with the CLI's display options: per-column
formatting by dtype (floats as `%.6f` with trailing zeros trimmed equally, or
`%.6e` when a value is below 1e-6 or the column is long and above 1e6; ints;
True/False; text with tab/CR/LF escaped; `NaN` for missing numbers and text,
`None` for Python's None in an object column), right-justified values and
headers, a space before numeric (int, float, bool) headers, columns joined by one
space, no truncation; `Empty DataFrame` / `Columns: [...]` / `Index: []` for no
rows. Widths count characters (Python's `len`).

**Comparison with Python** (`test/test_cli.jl`; fixture `test/data/python_cli/`
written by `uv run python julia/HTGEquilibrium/test/dump_cli.py`, 144 kB):

- `htg check` of the 5 case files, the 6 old case files, the template and 5
  template variants with problems (unknown species, aqueous phase above its
  water-density limit, 3500 °C, a mineral above its range, explicit points):
  **stdout, stderr and exit code identical**.
- `htg species` (16 variants: default, `--all`, `--phase` with plain text and
  regular expressions, `--search`, both, empty results, with the user species
  folder of `test/data/user_species`): **identical**; the five outputs above
  20 kB (`--all` variants, up to ~400 kB) are compared by sha256.
- `to_string` of 8 hand-made tables (every dtype, NaN, None, inf, -0.0, values
  needing exponent notation, text with line breaks, an empty table):
  **identical**, and the dtype inferred for each column is pandas'.
- 29 usage errors and help requests: **usage and error messages identical**
  (argparse 3.12 text: missing and unrecognized arguments, invalid choice,
  "expected one argument", "invalid int value", ignored explicit argument,
  abbreviations, `--`); the help texts below the usage differ on purpose (they
  say what the Julia version does not have yet).
- 9 errors reported as `Error: ...` (exit 2): missing files (Python's
  FileNotFoundError text, with the path as `pathlib` prints it), wrong data
  folder, duplicate case names, an invalid value, `new` into a missing folder;
  `new` on an existing file (exit 1): **identical**, except the TOML syntax
  error, whose text is Julia's TOML parser's (already accepted for case files).
- `htg new`: same message and the same file (the template with the name
  replaced; also `htg new sub.dir//other`).
- `htg run` of three independent points of `yakaboylu2013_fig7_8` (300, 389.6
  and 450 °C; the case adds 390 °C, its aqueous cutoff; 389.6 °C fails): exit
  code 1 as in Python; the output text before the summary (case line, preflight
  table, progress lines) identical; the summary table equal with numbers within
  1e-6 relative (the printed G has 7 digits); the output folder (file names,
  columns, rows, values, case.toml) equal by `test/compare_cases.jl`'s rules
  (all 4 points match). With `HTG_FULL_TESTS=1` the same run on 2 worker
  threads (with the plots note) gives the same files.

**Differences from Python** (deliberate):
- `htg plot` and `htg gui` are not available yet (M4): they parse their
  arguments as Python does, then print "Error: htg plot is not available in the
  Julia version yet ..." and exit with 2.
- `htg run` writes no plots: without `--no-plots` it says once, before the case
  lines, "Plots are not written by the Julia version yet (as with --no-plots);
  ...". Python printed a "Plot: ..." line per plot.
- `--workers N` means threads (above): "Starting N worker threads (the first
  results take about a minute while Julia compiles)" instead of "Starting N
  worker processes (each starts Julia; ...)", and the cap message.
- `--phase` with an invalid regular expression: "Error: --phase '[' is not a
  valid regular expression: ..." (exit 2); Python ended with a traceback (exit
  1). Patterns are PCRE (Julia), not Python's `re`; the same for ordinary
  patterns.
- An unexpected error (Python: a traceback, exit 1) prints Julia's error and
  backtrace, exit 1.
- `htg new` writes `\n` line ends on Windows too (Python's text mode wrote
  `\r\n`).
- Runs with a user species folder lack `user_species_used.toml` (see "Result
  tables ...": the run record comes with the user species port, M4).

**Manual check** (2026-11-24, from the shell, Linux): `julia/htg species
--search methane`, `julia/htg check cases/manure_gas_min.toml`, `julia/htg new
my_case.toml` (in an empty folder), `julia/htg plot results/x` (exit 2),
`julia/htg run` (usage error, exit 2) and
`julia/htg run cases/cellulose_matlab.toml cases/manure_gas_min.toml --out OUT --workers 4`
(24 points on 4 threads, all solved, exit 0; "Finished in 102 s" including
compilation, 147 s with Julia's start-up; tables printed as above).
`julia\htg.cmd` was not run (no Windows machine here). The Pkg app shim was
run once from a separate depot (`htg --data DIR species --phase ion`, and a
usage error with exit code 2).

## What remains for M2

Done: constants, formulas, water, HKF, band check; the data layer
(`databases.jl`, `saft.jl`, `registry.jl`, `user_species.jl` loading,
`standard_states.jl`, `model.jl` data part; `catalogue.json` and the
standard-state parts of `thermo.json`); case files and feed arithmetic
(`units.jl`, `case.jl` reading and writing, `feed_v2.jl`); the phase models
(`phases.jl`); the solver (`optimizer.jl`, `equilibrium.jl`, the data part of
`results.jl`, `Model.solve`; `solves.json` is tested); the result tables, CSV
files and choices (`tables.jl`, `results.jl`), sweeps and the preflight check
(`sweep.jl`), the worker pool on threads (`queue.jl`; BLAS on one thread); the
example cases run in both versions (`test/compare_cases.jl`: all 87 points
agree); the command line (`cli.jl`: `run`, `check`, `species`, `new`, the
launchers). `test/data/hkf_ions.toml` equals the registry's records (tested) and
can go when `test_water.jl` is switched over. Still to do:
- `htg plot` and `htg gui` (with the plotting and GUI ports, M4): today they
  exit with a message; `htg run` writes no plots.
- The user species run record of a saved run (`user_species_used.toml`, the
  note in case.toml; `user_species.run_record`), with the user species port (M4).
- Interim step (Python GUI calling the Julia core): the GUI uses `WorkerPool`
  (`submit`, `cancel`, `status`, `set_workers`, `restart`, `runs`, `Run`
  fields), `live_tables`/`load_saved` (`RunTables` with pandas DataFrames) and
  `Model.check`. Through juliacall: convert `DataTable`s to DataFrames (or let
  the GUI read the saved CSVs, which are identical), and do not call Python
  from the pool's worker threads (`on_point`/`on_done` run there): poll
  `status(pool)` and the `Run` fields from the GUI's timer instead, or pass
  callbacks that only push to a `Channel` that Python drains. Start Julia with
  threads (`PYTHON_JULIACALL_THREADS`, `PYTHON_JULIACALL_HANDLE_SIGNALS=yes`).
  Estimate (M2 CLI step): `runner.py`'s process pool (~300 of its 449 lines)
  becomes a juliacall wrapper of ~150 lines with the same API (cases passed
  as `case.to_toml()` text, equal in both versions; `Run` and `status()` as
  snapshots taken under the pool's lock; `wait` never on the event loop);
  `gui/app.py` (pool creation, start-up text) and `gui/runs.py` (`live_tables`
  from Julia tables or the saved CSVs) change by a few dozen lines; the GUI's
  data side (Python `Model`, species editing, preflight messages) stays
  Python; `src/htg/juliapkg.json` gets HTGEquilibrium as a dev path dependency;
  `tests/test_runner.py` is adapted. About 1–2 days, after a short spike that
  Julia threads keep solving (and GC works) while the Python thread sits in
  the GUI's event loop.
- Tests that do not depend on the old Python code (owner, PR #10): the sweep
  tests check element balances, signs, charge and carbon shares of every
  solution; still to add: KKT/stationarity and tangent-plane stability of
  absent phases as pass/fail checks, and comparisons with independent
  literature data, so that the Python references (`solves.json`, the Python
  dumps) can be retired.
- M5: cross-platform build of the Fortran library (JLL).

## Owner decisions (PR #4)

1. **Fortran source in git**, in its own folder, `vendor/supcrt92/`, with a
   README saying where and when it was copied.
2. **Band warning: warn only.** The defaults stay: "failed" means an implied
   density error above 1e-6, "near" means within 0.5 °C. The GUI does not move
   or drop temperatures.
3. **One water model throughout.** HKF takes every water property, including
   the 25 °C, 1 bar reference state, from the selected model; there is no
   option to mix models.
4. **Fortran `STOP`.** Two `STOP` statements in `power()` would end the whole
   app (never reached on valid input with the lock). A patched copy is
   acceptable; tracked for later in issue #5
   (https://github.com/ComKinBio/HTG-equilibrium/issues/5).
5. **Outside SUPCRT92's range** `water_properties` returns NaN, as pychnosz
   does. Fine for now.

## Owner decisions (PR #7)

1. **Gas product candidates:** fixed in Julia only (the Python app will be
   removed): NASA phase "G" or blank is a gas record.
2. **No bit matching:** plain Julia functions (`^`, `sum`); Python is a
   reference within tolerances (8e-14 is close).
3. **Water model per case:** a case key chosen with a button in the case
   editor; the model data hold none.

## Owner decisions (PR #9)

1. **Exact derivatives:** "No, just go straight to the exact derivatives."
   `hessian` is the ForwardDiff Jacobian of mu for fluid and electrolyte
   phases, the only path for the optimizer; mu from the fugacity coefficients
   stays the gradient. Python's forward-difference Hessian is not ported.
2. **Water by name:** "Yes, make it check water specifically so it is not
   dependent on the ordering of species." Every emptiness check of the
   electrolyte phase uses the water amount (tested with water not listed first,
   and with water not selected).

## Owner decisions (PR #10)

1. **Dropped solvers:** "Yes, drop them for good." `trust-constr` and
   `slsqp-fd` are removed (with `fd_ftol`, `fd_eps`); old cases that name them
   load as `ipopt` with a note in `case.legacy_notes` and are written as
   `ipopt` (see "Case files: reading").
2. **Reference points of `slsqp-fd`:** "Yes, the stored reference used a bad
   optimizer. So those should be treated as a rough guide. In the future, we
   should think about making testing a bit more independent of the old python
   code." They are solved with IPOPT; a lower G is accepted when it matches the
   Python app's IPOPT solution (`LOWER_THAN_FD`). Independent tests are listed
   under "What remains for M2".
3. **HiGHS_jll 1.12.0:** "Yes, that's fine." Pinned to SciPy's HiGHS so the
   LP starts equal Python's.
4. **SLSQP limit:** "Multiply the default by 10." SLSQP may use
   10 × `max_iterations` function evaluations (50000 by default; see "Solver",
   differences).
