# HTGEquilibrium.jl

The Julia science core of the HTG equilibrium app (Phase 2, milestone M2 of
`docs/spec/julia_migration_plan.md`). It is ported from the Python package
`src/htg` one module at a time, keeping the Python names so the two are easy to
compare. So far: constants, formula parsing, water properties, HKF standard
states of aqueous ions, the SUPCRT92 failure-band check that the Julia GUI
will use to warn about sweeps, the data layer (databases, species
registry, user species folder, standard states, model loading) and reading
case files (units, feed arithmetic).

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
| `test/dump_cases.py`, `test/python_cases.jl` | Write / compare `test/data/cases_python.json` (~0.5 MB, not committed) | – |
| `test/data_setup.jl` | Shared setup of the three tests below: the reference models, Python's `same`/`close` | `tests/conftest.py`, `tests/helpers.py` |
| `test/test_databases.jl` | Readers and Python helpers (CSV dtypes, XML, NASA-7, ordered TOML) | – |
| `test/test_registry.jl` | `catalogue.json` (all of it), a user species folder, shipped `data/` smoke test, threaded reads | `tests/test_reference_catalogue.py`, `tests/test_user_species.py` (loading) |
| `test/test_standard_states.jl` | `thermo.json`: molecular and solid g°, SAFT parameters, mu0 of the fluid and electrolyte phases | `tests/test_reference_thermo.py` |
| `test/data/user_species/`, `test/data/user_species_python.json` | A user species folder written by the Python app and what Python loads from it; written by `test/dump_user_species.py` | – |

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

**Comparison with Python.** `uv run python julia/HTGEquilibrium/test/dump_cases.py`
writes `test/data/cases_python.json` (~0.5 MB, ignored by git; without it the
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
such key; every Python case has the default).

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
  editing/importing user species (M4), `Model.solve/sweep/check`.

**For the phases step.** `resolve_parameters` gives the parameter sets
(`m`, `sigma`, `epsilon`, `epsilonAB`, `kappaAB`, `association_scheme`, `Mw`)
that `build_pcsaft_model` needs; components are `all_molecular(registry)[name]["formula"]`
with `registry.binary_interactions`. The electrolyte phase needs
`aqueous_reference_neutrals/ions`, `probe_weights`, `ion_charge`,
`registry.ion` (formula, elements), `ion_epcsaft_table` (m, sigma, epsilon,
dielectric_constant), `hkf_ion_g0(water_model, model, ions, T, P)`,
`nasa_g0(model, neutrals, T)` and `water_properties(water_model, T, P)` for the
validity check and the dielectric constant (Python calls `supcrt_water` there;
Julia takes the case's water model). Still untested from `thermo.json`: the fluid-phase G and ln phi
points (rel 1e-10 / abs 1e-8 and rtol 1e-9 / atol 1e-11), and the electrolyte
`lnphi_infinity`, G and ln phi.

## What remains for M2

Done: constants, formulas, water, HKF, band check; the data layer
(`databases.jl`, `saft.jl`, `registry.jl`, `user_species.jl` loading,
`standard_states.jl`, `model.jl` data part; `catalogue.json` and the
standard-state parts of `thermo.json`); case files and feed arithmetic
(`units.jl`, `case.jl` reading and writing, `feed_v2.jl`). `test/data/hkf_ions.toml`
equals the registry's records (tested) and can go when `test_water.jl` is
switched over. Still to do:
- `phases.jl` (electrolyte phase with the water model selectable; validity
  check `rho >= 350 and epsilon > 1`), `equilibrium.jl`, `optimizer.jl`,
  `results.jl`, `sweep.jl`, `queue.jl`, `cli.jl`; the fluid-phase and
  electrolyte G / ln phi parts of `thermo.json`; `solves.json` tests; running
  the example cases in both versions.
- Passing the case's water model (`water_model_of(case)`, see "Case files")
  to phases; results record it.
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
