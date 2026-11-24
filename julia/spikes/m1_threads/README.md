# M1 spike: thread safety and parallelism of the Julia solver

Milestone M1 of `docs/spec/julia_migration_plan.md`: can one Julia process solve
several equilibrium points at once on threads, calling Clapeyron and Ipopt the way
the app does, with the same results every time?

**Outcome: Threads.** One process, a task queue, N solver threads. No thread-safety
problem was found in 8 640 solves on 4–8 threads, 3 240 of them sharing one set
of Clapeyron models (every result bitwise identical to the serial run), and 8 threads use about 2.2 GB in total against 11.5 GB for 8
Distributed workers. Distributed workers are about 30 % faster at 8-way
parallelism today because each worker has its own garbage collector; M2 should
cut allocations in the objective and Hessian to close that gap (see
"Problems found").

## Method

### What was replayed

This is a **real replay**, not a representative problem. For every solved case in
`tests/reference/data/solves.json` (21 points; run with the app's default solver,
IPOPT; 27 (point, start) solves including the 7-start multistart case),
`dump_problems.py` captures from the running Python app everything that crosses
into Julia for one point:

- the Clapeyron model inputs: `PCSAFT(components; userlocations)` per fluid phase,
  and for the electrolyte phase `ePCSAFT(components, charges, ideal, PCSAFT,
  hsdDH(components, ConstRSP(eps_r), refs), refs)`, with exactly the
  parameters the app builds (recorded by wrapping `pyclapeyron.PCSAFT`);
- the phase data Python computes (NASA/HKF standard states, T, P, root, probe
  weights, positions) and the IPOPT problem exactly as `optimizer._ipopt` builds
  it (A, b, inequalities, bounds, Hessian structure, options, starting points);
- check values computed by Python through juliacall (fugacity coefficients,
  objective, gradient and Hessian at the start, ln φ at infinite dilution), and
  the Python app's raw IPOPT result from every start.

`replay.jl` ports the Julia-facing code line for line: `FluidPhase`,
`ElectrolytePhase` (electroneutral probes, infinite-dilution extrapolation),
`_hessian` (forward differences of ln φ), `GibbsProblem.fun_and_grad` /
`hessian_blocks`, and `_ipopt` + `_JULIA_IPOPT` (same callbacks, options,
NaN-on-EOS-failure). The only Clapeyron call the app makes while solving is
`fugacity_coefficient(model, P, T, n; phase = "vapour" | "liquid")`; the replay
makes the same call. Each task builds its own Clapeyron models and its own Ipopt
problem, as the app does per point.

Not replayed (pure Python/SciPy today, ported in M2): LP starts (taken from the
dump), the SLSQP fallback and polish, KKT and tangent-plane diagnostics, and the
standard-state calculations (pychnosz/HKF; values taken from the dump). Their
Julia ports (HiGHS.jl, NLopt.jl, own water/HKF code) need the same one-object-per-
task rule; see "For M2".

### Runs

| script | what |
|---|---|
| `dump_problems.py` | `uv run python julia/spikes/m1_threads/dump_problems.py` → `problems.json` |
| `run_threads.jl` | `julia --project=. -t N run_threads.jl` (N = 1, 4, 8) |
| `run_distributed.jl` | `M1_WORKERS=N julia --project=. run_distributed.jl` (N = 4, 8) |
| `check_phi_bits.jl` | explains the last-bit differences from Python (below) |

`run_threads.jl`: warm-up pass (compilation); checks against the Python values;
3 serial runs (baseline); then **20 runs of a batch of 108 solves** (4 copies of
the 27 tasks, spawned at once with `Threads.@spawn` in a new random order each
run); then a **shared-model stress test** (5 rounds; for each of the 27 tasks,
one set of Clapeyron models used by 2N concurrent solves). Every solve is
compared bitwise (x, objective, status, Hessian-call count, function-call count)
with the serial baseline, and to the reference tolerances with Python's result.
`run_distributed.jl` runs the same 20 × 108 batches with `pmap` on N workers
(1 thread each), the master holding no solver, and compares bitwise with the
serial baseline of `threads_1.json`. BLAS is set to 1 thread everywhere (see
"Problems found"). Raw output is in `results/`.

### Versions

Julia 1.12.7; the spike's `Project.toml`/`Manifest.toml` pin the same versions as
`src/htg/juliapkg.json` (manifest taken from the app's Julia environment, minus
PythonCall): Clapeyron 0.6.30, Ipopt.jl 1.16.0, Ipopt_jll 300.1400.1902 =
**Ipopt 3.14.19**, MUMPS_seq_jll 500.900.100 = **MUMPS 5.9.1** (sequential), linear
solver **MUMPS** (Ipopt default), OpenBLAS32_jll 0.3.34, SPRAL_jll 2025.9.18
(present, not used). Machine: Intel i7-1370P (6 performance + 8 efficiency cores,
20 threads), 29 GB, Linux; other jobs were running at the same time (load average
2–9), so timings are indicative.

## Results

### Replay is faithful to the Python app

| check (21 points) | result |
|---|---|
| fugacity-coefficient vectors bitwise equal to Python/juliacall | 48 / 53; the other 5 differ by ≤ 5.4e-14 relative |
| ln φ∞ of ions | max diff 1.3e-14 |
| objective at start | max diff 1.9e-9 J/mol C |
| gradient at start | max diff 1.1e-16 (scaled) |
| Hessian blocks at start | max rel diff 2.6e-12 |
| IPOPT result vs Python IPOPT from the same start (G rel 1e-9, amounts rtol 1e-5 / atol 1e-7, same status) | **27 / 27 pass**; worst ΔG 3.7e-9 J/mol C, worst Δx 3.9e-13 mol |
| same number of IPOPT Hessian calls as Python | 22 / 27 (others ±1–2 iterations, same answer) |

The last-bit differences come from the array type, not from Clapeyron:
juliacall passes numpy arrays as `PyArray` (not a dense Julia array), so
Clapeyron's sums take Julia's generic path; with a plain `AbstractVector`
wrapper the replay matches Python bitwise for 51 / 53 (`check_phi_bits.jl`; the
last 2 differ in the probe composition, numpy vs Julia `dot`). Against the stored
reference data: the 7 IPOPT-reference points match to ≤ 4e-16 relative, except
`full_text_250C_ipopt`, which reaches the lower `full_text_250C` solution as the
Python test expects (`IMPROVED`); for the 14 points whose reference comes from the
finite-difference SLSQP, IPOPT (Python and Julia alike) reaches the same or a lower
G (`cellulose_400C` 2.2e-6 relative lower, the others within 7e-9).

### Thread safety

| run | solves | differ bitwise from serial | outside tolerance | errors |
|---|---|---|---|---|
| serial, 3 repeats | 81 | 0 | 0 | 0 |
| 1 thread, 20 × 108 | 2 160 | 0 | 0 | 0 |
| 4 threads, 20 × 108 | 2 160 | 0 | 0 | 0 |
| 8 threads, 20 × 108 | 2 160 | 0 | 0 | 0 |
| 8 threads (`--gcthreads=8,1`), 10 × 108 | 1 080 | 0 | 0 | 0 |
| shared Clapeyron models, 4 threads (8 concurrent per model) | 1 080 | 0 | – | 0 |
| shared Clapeyron models, 8 threads (16 concurrent per model) | 2 160 | 0 | – | 0 |
| Distributed, 4 workers, 20 × 108 | 2 160 | 0 | – | 0 |
| Distributed, 8 workers, 20 × 108 | 2 160 | 0 | – | 0 |

**Pass.** Results are deterministic: serial, threaded and Distributed solves are
bitwise identical, including iteration counts, in every run.

What was checked in the libraries:

- **Ipopt** (one `IpoptProblem` per task): Ipopt.jl keeps no global state; the C
  callbacks get the problem through `user_data` and run on the calling thread.
  Ipopt 3.14.19 serialises every call into MUMPS with a process-wide mutex
  (`mumps_call_mutex` in `libipopt.so.3.14.19`, `IpMumpsSolverInterface`), because
  sequential MUMPS is not thread-safe. So concurrent solves are safe; only the
  factorisations take turns, which costs little here (≤ 47 variables).
- **Clapeyron** models are immutable parameter structs; `fugacity_coefficient`
  with a given phase uses one volume root and no caches. The only internal
  threading (`StableTasks.@spawn` in `default_volume_impl`) runs only for
  `phase = :unknown`, which the app never uses. The only mutable global
  (`PARAM_LOCATION`, scratch space for `ParamTable`) is not touched by model
  construction from `userlocations`. Model construction on 8 threads at once and
  16 solves sharing one model both gave bitwise identical results.
- **Our code**: the only cache is `FluidPhase._pure_mu` (per phase object, so
  per point/task), plus the IPOPT `evaluate` last-x cache (per solve). The Python
  module-level `_ipopt_solve` (a compiled function) has no Julia counterpart. Rule
  for M2: one phase object and one Ipopt problem per task; Clapeyron models may be
  shared but there is no need to.

### Speed

Wall time for one batch of 108 solves (median of 20 runs; serial batch in a
1-thread process: 2.94 s):

| configuration | median s / batch | speed-up | GC share of wall |
|---|---|---|---|
| serial (1 thread) | 2.94 | 1.0× | 12 % |
| Threads, 4 | 1.29 | 2.3× | 21 % |
| Threads, 8 | 1.18 | 2.5× | 24 % |
| Threads, 8, `--gcthreads=8,1` (10 runs) | 1.17 | 2.5× | 20 % |
| Distributed, 4 workers | 1.20 | 2.5× | – |
| Distributed, 8 workers | 0.91 | 3.2× | – |

Serially, the Julia replay takes 0.77 s for the 26 tasks after the first, against
1.37 s for the Python app's IPOPT through juliacall from the same starts.

Start-up without a sysimage: one process loads in a few seconds and compiles in
the first pass (63–73 s); Distributed workers load in 9–12 s but compile in each
worker (125 s for 4, 200 s for 8 in parallel, competing for cores).

### Memory

After the same work (warm-up, checks, 20 batches), measured from `/proc`:

| configuration | processes | RSS total | PSS total | peak RSS (sum) |
|---|---|---|---|---|
| Threads, 1 | 1 | 1.84 GB | 1.83 GB | 1.91 GB |
| Threads, 4 | 1 | 2.08 GB | 1.99 GB | 2.09 GB |
| Threads, 8 | 1 | 2.16 GB | 2.06 GB | 2.36 GB |
| Distributed, 4 workers | 1 + 4 | 5.99 GB (workers 5.50) | 5.09 GB | 6.50 GB |
| Distributed, 8 workers | 1 + 8 | 11.54 GB (workers 11.04) | 9.76 GB | 12.44 GB |

Each solver process is about 1.35–1.4 GB (Julia + Clapeyron + Ipopt + JIT-compiled
code; live Julia heap only ~150 MB, the rest is code and heap the GC keeps). Extra
threads cost about 0.1 GB each; extra workers cost the full 1.4 GB each, close
to today's 1.5 GB Python + Julia worker. The master of the Distributed runs (no
solver loaded) is 0.5 GB.

## Problems found

1. **BLAS threads.** By default OpenBLAS starts 10 threads, which Ipopt/MUMPS use
   through libblastrampoline. For these small problems that made the serial run
   slower (1.19 s vs 0.86 s for the 27 tasks in a trial) and oversubscribes cores
   when solver threads run. The app should call `BLAS.set_num_threads(1)`.
2. **GC limits thread speed-up.** A solve allocates ~30 MB (3.3 GB per batch of
   108), mostly in the finite-difference Hessian (n + 1 fugacity calls per phase,
   each allocating) and the vectorised objective. Julia's GC stops all threads,
   so its share grows from 12 % (serial) to 24 % (8 threads), and 8 threads give
   2.5× where 8 independent workers give 3.2×. More GC threads did not help. M2
   should reduce allocations (preallocated buffers, `fugacity_coefficient!`, or an
   analytic/ForwardDiff Hessian) before tuning further. On this hybrid CPU only 6
   cores are performance cores, which also caps scaling at 8.
3. **MUMPS is serialised by Ipopt** (mutex above). Not a correctness issue; it
   would matter only for much larger problems. SPRAL (shipped in Ipopt_jll) is
   the alternative if that ever becomes a bottleneck.
4. **Not bitwise equal to the Python app** in 5 of 53 fugacity vectors and in the
   iteration count of 5 of 27 solves (array type and summation order, see above).
   All results agree with Python to the reference tolerances, but the Julia port
   cannot be validated by bitwise comparison with Python; use the tolerances of
   `tests/test_reference_solves.py`.

## Decision for the owner

- **Recommended: Threads.** One Julia process with a task queue and N solver
  threads (N ≤ performance cores), BLAS at 1 thread, the GUI on its own thread.
  It passed every thread-safety check and uses ~2.2 GB for 8 threads instead of
  ~11.5 GB (Distributed) or ~12 GB (8 Python workers today).
- Accept that, until M2 reduces allocations, threads are somewhat slower than
  the same number of Distributed workers (here 1.18 s vs 0.91 s per 108 solves
  at 8-way). Distributed stays a fallback that needs no code change in the solver
  (the same `solve_task` runs in both), at 1.4 GB per worker.

## For M2

- Keep the replay structure: a pure function per point (build phases, build
  Ipopt problem, solve) with no global mutable state.
- New Julia dependencies (HiGHS.jl for LP starts, NLopt.jl for SLSQP, the
  water/HKF code) must also be created per task; rerun `run_threads.jl` with the
  full solve once they exist.
- Ship a sysimage: it removes the 60–70 s first-pass compilation and should also
  lower the per-process code memory measured above.

## Files

`Project.toml`, `Manifest.toml` (pinned environment), `dump_problems.py`,
`problems.json` (captured inputs, 280 kB), `replay.jl`, `memory.jl`,
`run_threads.jl`, `run_distributed.jl`, `check_phi_bits.jl`, `results/`
(JSON summaries, logs, `/usr/bin/time -v` output).
