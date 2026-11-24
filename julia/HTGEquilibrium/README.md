# HTGEquilibrium.jl

The Julia science core of the HTG equilibrium app (Phase 2, milestone M2 of
`docs/spec/julia_migration_plan.md`). It is ported from the Python package
`src/htg` one module at a time, keeping the Python names so the two are easy to
compare. So far: constants, formula parsing, water properties, HKF standard
states of aqueous ions, and the SUPCRT92 failure-band check that the Julia GUI
will use to warn about sweeps.

| File | Contents | Python counterpart |
|---|---|---|
| `src/constants.jl` | R, P0, reference states, atomic weights, phase names; `pysum` (Python's `sum()` bit for bit) | `constants.py` |
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

**Matching Python's arithmetic.** Python 3.12's `sum()` of floats uses
compensated summation; Julia's `sum` does not, so results can differ in the last
bit. Use `pysum` where the Python code uses `sum()` and the result is compared
with Python. Python's `dict` keeps insertion order; use `OrderedDict` where order
matters (element lists, species lists).

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

## What remains for M2

This package does not yet have the rest of `standard_states.py` or anything
that uses it:
- `constants.jl`, `databases.jl` (NASA-7 parser, CSV SAFT tables, EzXML for
  `supcrt92.xml`) and `registry.jl`. Then the HKF records come from
  `data/supcrt92.xml` and `test/data/hkf_ions.toml` can go.
- `standard_states.jl`: `nasa_g0`, `supcrt_mineral_properties`, `solid_g0`, and
  the remaining `thermo.json` tests (molecular, solid, SAFT parameters).
- `phases.jl` (electrolyte phase with the water model selectable; validity
  check `rho >= 350 and epsilon > 1`), `equilibrium.jl`, `optimizer.jl`,
  `case.jl`, `results.jl`, `sweep.jl`, `queue.jl`, `cli.jl`; `solves.json` and
  `catalogue.json` tests; running the example cases in both versions.
- Where the user chooses the water model (case file key, GUI setting) and how
  results record it.
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
