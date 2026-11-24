# M1 thread-safety test: serial vs Threads.@spawn, repeated.
#
#   julia --project=. -t N run_threads.jl
#   environment: M1_REPS (threaded runs, default 20), M1_MULT (copies of the task
#   list per run, default 4), M1_BLAS (BLAS threads, default 1), M1_OUT (json path)
#
# 1. checks the replay against the values Python computed through juliacall;
# 2. serial baseline (3 runs, must be bitwise identical to each other) and
#    comparison with the Python app's IPOPT result from the same start;
# 3. M1_REPS threaded runs: M1_MULT copies of every (point, start) task of
#    problems.json spawned at once in shuffled order; each task builds its own
#    Clapeyron models and Ipopt problem. Every result is compared bitwise with
#    the serial baseline and to the reference tolerances with Python's result;
# 4. shared-model stress: per point, the Clapeyron models are built once and
#    2*nthreads tasks solve that point concurrently with the same model objects;
# 5. memory (RSS/PSS/peak) of this process.

using Printf, Random, Statistics, LinearAlgebra, JSON
include(joinpath(@__DIR__, "replay.jl"))
using .Replay
using Ipopt
include(joinpath(@__DIR__, "memory.jl"))

const REPS = parse(Int, get(ENV, "M1_REPS", "20"))
const MULT = parse(Int, get(ENV, "M1_MULT", "4"))
const NBLAS = parse(Int, get(ENV, "M1_BLAS", "1"))
const OUT = get(ENV, "M1_OUT", joinpath(@__DIR__, "results", "threads_$(Threads.nthreads()).json"))

# ---------------------------------------------------------------- replay vs Python check values
function check_against_python(cases)
    worst = Dict("phi_rel" => 0.0, "lnphi_inf" => 0.0, "f0_J" => 0.0, "grad0" => 0.0, "hess0_rel" => 0.0)
    nphi, nphi_bitwise = 0, 0
    for c in cases
        phases = Replay.build_phases(c)
        for (name, ph) in phases
            rec = c["phases"][name]
            if haskey(rec, "check_n")
                n = Float64.(rec["check_n"])
                φ = ph isa Replay.FluidPhase ? Replay.phis(ph, n) : Replay.phi(ph, n)
                ref = Float64.(rec["check_phi"])
                nphi += 1
                nphi_bitwise += φ == ref
                worst["phi_rel"] = max(worst["phi_rel"], maximum(abs.(φ ./ ref .- 1)))
            end
            if ph isa Replay.ElectrolytePhase
                worst["lnphi_inf"] = max(worst["lnphi_inf"],
                                         maximum(abs.(ph.lnphi_infinity .- Float64.(rec["lnphi_infinity"]))))
            end
        end
        pr = Replay.GibbsProblem(c, phases)
        x0 = Float64.(c["x_start"])
        f, g = Replay.fun_and_grad(pr, x0)
        worst["f0_J"] = max(worst["f0_J"], abs(f - c["check"]["f0"]) * c["scale"])
        worst["grad0"] = max(worst["grad0"], maximum(abs.(g .- Float64.(c["check"]["grad0"]))))
        names = [n for (n, p) in phases if !(p isa Replay.SolidPhase)]
        for (k, (pos, H)) in enumerate(Replay.hessian_blocks(pr, x0))
            ref = c["check"]["hessian_blocks0"][names[k]]
            @assert pos == Int.(ref["pos"]) .+ 1
            Href = reduce(vcat, permutedims(Float64.(r)) for r in ref["H"])
            worst["hess0_rel"] = max(worst["hess0_rel"], maximum(abs.(H .- Href)) / max(1.0, maximum(abs.(Href))))
        end
    end
    return nphi, nphi_bitwise, worst
end

# ---------------------------------------------------------------- comparisons
G_of(c, r) = c["G0"] + r.obj * c["scale"]

"""Julia result vs the Python app's raw IPOPT result from the same start, with the
tolerances of tests/test_reference_solves.py (G rel 1e-9; amounts rtol 1e-5, atol 1e-7)."""
function vs_python(c, s, r)
    py = c["python_ipopt"][s]
    Gj, Gp = G_of(c, r), c["G0"] + py["obj"] * c["scale"]
    xp = Float64.(py["x"])
    ok = abs(Gj - Gp) <= 1e-9 * abs(Gp) && all(abs.(r.x .- xp) .<= 1e-7 .+ 1e-5 .* abs.(xp)) &&
         r.status == py["status"]
    return (ok = ok, dG = abs(Gj - Gp), dx = maximum(abs.(r.x .- xp)), same_hess = r.hessian_evals == py["hessian_evals"])
end

same(a, b) = a.x == b.x && a.obj === b.obj && a.status == b.status && a.hessian_evals == b.hessian_evals &&
             a.nfev == b.nfev && a.eos_failures == b.eos_failures

run_serial(cases, tl) = (t0 = time(); r = [Replay.solve_task(cases[i], s) for (i, s) in tl]; (r, time() - t0))

function run_threaded(cases, tl, rng)
    order = shuffle(rng, collect(eachindex(tl)))
    t0 = time()
    futures = Vector{Task}(undef, length(tl))
    for k in order
        futures[k] = Threads.@spawn Replay.solve_task(cases[tl[k][1]], tl[k][2])
    end
    res = [fetch(f) for f in futures]
    return res, time() - t0
end

function main()
    BLAS.set_num_threads(NBLAS)
    cases = Replay.load_cases(joinpath(@__DIR__, "problems.json"))
    tasklist = Replay.tasks(cases)
    batch = repeat(tasklist, MULT)
    println("Julia $(VERSION), threads = $(Threads.nthreads()), BLAS threads = $(BLAS.get_num_threads()), ",
            "Ipopt $(Ipopt.GetIpoptVersion()); $(length(cases)) points, $(length(tasklist)) (point, start) tasks; ",
            "batch = $MULT x $(length(tasklist)) = $(length(batch)) solves")
    S = Dict{String,Any}("threads" => Threads.nthreads(), "blas_threads" => NBLAS, "reps" => REPS, "mult" => MULT,
                         "julia" => string(VERSION), "ipopt" => Ipopt.GetIpoptVersion(),
                         "ntasks" => length(tasklist), "batch" => length(batch))

    println("\n[0] warm-up (first serial pass, includes compilation)")
    S["memory_start"] = memory_report()
    twarm = @elapsed run_serial(cases, tasklist)
    @printf("    %.1f s\n", twarm)
    S["warmup_seconds"] = twarm

    println("[1] replay vs Python check values (same Clapeyron inputs, Python results via juliacall)")
    nphi, nbit, worst = check_against_python(cases)
    println("    fugacity-coefficient vectors bitwise equal to Python's: $nbit / $nphi")
    for (k, v) in sort(collect(worst)); @printf("    max |diff| %-10s %.3e\n", k, v); end
    S["check_phi_bitwise"] = [nbit, nphi]
    S["check_worst"] = worst

    println("[2] serial baseline")
    base, t1 = run_serial(cases, tasklist)
    serial_task_times = [t1]
    repeatable = true
    for _ in 1:2
        r, t = run_serial(cases, tasklist)
        push!(serial_task_times, t)
        repeatable &= all(same.(r, base))
    end
    # serial time of one batch (same work as one threaded run)
    gc0 = Base.gc_num()
    _, tbatch = run_serial(cases, batch)
    gcs = Base.GC_Diff(Base.gc_num(), gc0)
    @printf("    serial, %d tasks: %s s; repeats bitwise identical: %s\n", length(tasklist),
            join([@sprintf("%.2f", t) for t in serial_task_times], ", "), repeatable)
    @printf("    serial, one batch of %d solves: %.2f s (GC %.0f%%, %.0f MB allocated)\n", length(batch), tbatch,
            100 * gcs.total_time / 1e9 / tbatch, gcs.allocd / 2^20)
    vp = [vs_python(cases[i], s, base[k]) for (k, (i, s)) in enumerate(tasklist)]
    println("    vs Python IPOPT from the same start (G rel 1e-9, x rtol 1e-5/atol 1e-7, same status): ",
            count(v -> v.ok, vp), "/", length(vp), " pass; same IPOPT Hessian-call count: ",
            count(v -> v.same_hess, vp), "/", length(vp))
    @printf("    worst |dG| = %.2e J/mol C, worst |dx| = %.2e mol\n", maximum(v -> v.dG, vp), maximum(v -> v.dx, vp))
    S["serial_seconds_tasks"] = serial_task_times
    S["serial_seconds_batch"] = tbatch
    S["serial_gc_fraction"] = gcs.total_time / 1e9 / tbatch
    S["serial_repeatable"] = repeatable
    S["pass_vs_python"] = [count(v -> v.ok, vp), length(vp)]
    S["per_task"] = [Dict("case" => cases[i]["name"], "start" => cases[i]["python_ipopt"][s]["start"],
                          "n" => length(base[k].x), "status" => base[k].status, "G_julia" => G_of(cases[i], base[k]),
                          "G_python" => cases[i]["G0"] + cases[i]["python_ipopt"][s]["obj"] * cases[i]["scale"],
                          "G_python_full_solve" => cases[i]["python_solve"]["G"],
                          "G_reference" => cases[i]["reference"]["G"],
                          "reference_solver" => cases[i]["reference_solver"],
                          "dG_python" => vp[k].dG, "dx_python" => vp[k].dx, "pass_python" => vp[k].ok,
                          "hess_julia" => base[k].hessian_evals,
                          "hess_python" => cases[i]["python_ipopt"][s]["hessian_evals"],
                          "eos_failures" => base[k].eos_failures, "seconds_julia" => base[k].seconds,
                          "seconds_python_ipopt" => cases[i]["python_ipopt"][s]["seconds"],
                          "x" => base[k].x, "obj" => base[k].obj, "nfev" => base[k].nfev)
                     for (k, (i, s)) in enumerate(tasklist)]

    println("[3] threaded: $REPS runs x $(length(batch)) solves on $(Threads.nthreads()) threads")
    rng = MersenneTwister(1)
    times, nbit_fail, ntol_fail, nerr = Float64[], 0, 0, 0
    basebatch = repeat(base, MULT)
    gc0 = Base.gc_num()
    for rep in 1:REPS
        res = try
            r, t = run_threaded(cases, batch, rng)
            push!(times, t)
            r
        catch e
            nerr += 1
            println("    run $rep: ERROR ", first(sprint(showerror, e), 400))
            continue
        end
        bad = findall(.!same.(res, basebatch))
        tol = [k for k in eachindex(res) if !vs_python(cases[batch[k][1]], batch[k][2], res[k]).ok]
        nbit_fail += length(bad)
        ntol_fail += length(tol)
        isempty(bad) && isempty(tol) ||
            println("    run $rep: $(length(bad)) differ bitwise, $(length(tol)) outside tolerance: ",
                    unique([cases[batch[k][1]]["name"] for k in bad]))
    end
    gct = Base.GC_Diff(Base.gc_num(), gc0)
    @printf("    wall time per batch: median %.2f s (min %.2f, max %.2f); speed-up vs serial batch %.2fx; GC %.0f%% of wall\n",
            median(times), minimum(times), maximum(times), tbatch / median(times),
            100 * gct.total_time / 1e9 / sum(times))
    println("    solves differing bitwise from serial: $nbit_fail / $(REPS * length(batch)); ",
            "outside tolerance: $ntol_fail; errors: $nerr")
    S["threaded_seconds"] = times
    S["threaded_speedup_median"] = tbatch / median(times)
    S["threaded_gc_fraction"] = gct.total_time / 1e9 / sum(times)
    S["threaded_solves"] = REPS * length(batch)
    S["threaded_bitwise_failures"] = nbit_fail
    S["threaded_tolerance_failures"] = ntol_fail
    S["threaded_errors"] = nerr

    nshared = 2 * Threads.nthreads()
    shared_reps = max(1, REPS ÷ 4)
    println("[4] shared-model stress: $shared_reps x $nshared concurrent solves of each task with one set of Clapeyron models")
    sfail, sruns = 0, 0
    for _ in 1:shared_reps, (k, (i, s)) in enumerate(tasklist)
        models = Replay.build_models(cases[i])
        fs = [Threads.@spawn Replay.solve_task(cases[i], s; shared = models) for _ in 1:nshared]
        for f in fs
            sruns += 1
            r = try fetch(f) catch e; println("    ERROR $(cases[i]["name"]): ", first(sprint(showerror, e), 300)); nothing end
            (r === nothing || !same(r, base[k])) && (sfail += 1)
        end
    end
    println("    shared-model solves differing bitwise from serial: $sfail / $sruns")
    S["shared_model_runs"] = sruns
    S["shared_model_failures"] = sfail

    GC.gc(true)
    S["memory"] = memory_report()
    println("[5] memory of this process: ", S["memory"])
    mkpath(dirname(OUT))
    open(io -> JSON.json(io, S; pretty = true), OUT, "w")
    println("wrote $OUT")
end

main()
