# M1 comparison: the same batches as run_threads.jl on N Distributed workers.
#
#   julia --project=. run_distributed.jl
#   environment: M1_WORKERS (default 4), M1_REPS (default 20), M1_MULT (default 4),
#   M1_BASELINE (a results/threads_*.json whose serial results are the reference),
#   M1_OUT (json path)
#
# The master process only hands out tasks (like a GUI process would): it does
# not load Clapeyron or Ipopt. Each worker (one thread, BLAS 1 thread) loads
# replay.jl and the problems. Results are compared bitwise with the serial
# baseline of a threads run; memory is summed over master + workers.

using Distributed, Printf, Statistics, Random, JSON
include(joinpath(@__DIR__, "memory.jl"))

const N = parse(Int, get(ENV, "M1_WORKERS", "4"))
const REPS = parse(Int, get(ENV, "M1_REPS", "20"))
const MULT = parse(Int, get(ENV, "M1_MULT", "4"))
const BASELINE = get(ENV, "M1_BASELINE", joinpath(@__DIR__, "results", "threads_1.json"))
const OUT = get(ENV, "M1_OUT", joinpath(@__DIR__, "results", "distributed_$(N).json"))

const T0 = time()
addprocs(N; exeflags = ["--project=$(@__DIR__)", "--threads=1"])
const DIR = @__DIR__
@everywhere workers() begin
    using LinearAlgebra
    BLAS.set_num_threads(1)
    include(joinpath($DIR, "replay.jl"))
    include(joinpath($DIR, "memory.jl"))
    const CASES = Replay.load_cases(joinpath($DIR, "problems.json"))
    solve_k(i, s) = (r = Replay.solve_task(CASES[i], s); (r.x, r.obj, r.status, r.hessian_evals, r.nfev))
    warm() = (for (i, s) in Replay.tasks(CASES); Replay.solve_task(CASES[i], s); end)
end
const START_SECONDS = time() - T0

function main()
    S = Dict{String,Any}("workers" => N, "reps" => REPS, "mult" => MULT)
    dir = DIR
    S["start_seconds"] = START_SECONDS
    t1 = time()
    @sync for w in workers()
        @async remotecall_fetch(() -> Main.warm(), w)
    end
    S["warmup_seconds"] = time() - t1
    @printf("%d workers: started and loaded in %.1f s, warm-up (compilation, in parallel) %.1f s\n",
            N, S["start_seconds"], S["warmup_seconds"])

    base = JSON.parsefile(BASELINE)["per_task"]
    cases = JSON.parsefile(joinpath(dir, "problems.json"))["cases"]
    tasklist = [(i, s) for (i, c) in enumerate(cases) for s in eachindex(c["python_ipopt"])]
    @assert length(base) == length(tasklist)
    batch = repeat(collect(eachindex(tasklist)), MULT)
    rng = MersenneTwister(1)
    times, nfail = Float64[], 0
    for rep in 1:REPS
        order = shuffle(rng, collect(eachindex(batch)))
        t = @elapsed res = pmap(k -> Main.solve_k(tasklist[batch[k]]...), order)
        push!(times, t)
        for (k, r) in zip(order, res)
            b = base[batch[k]]
            same = r[1] == Float64.(b["x"]) && r[2] === Float64(b["obj"]) && r[3] == b["status"] &&
                   r[4] == b["hess_julia"] && r[5] == b["nfev"]
            nfail += !same
        end
    end
    @printf("wall time per batch of %d solves: median %.2f s (min %.2f, max %.2f)\n",
            length(batch), median(times), minimum(times), maximum(times))
    println("solves differing bitwise from the serial baseline: $nfail / $(REPS * length(batch))")
    S["seconds"] = times
    S["bitwise_failures"] = nfail
    S["solves"] = REPS * length(batch)

    for w in workers(); remotecall_fetch(() -> GC.gc(true), w); end
    master = memory_report()
    ws = [remotecall_fetch(() -> Main.memory_report(), w) for w in workers()]
    total(key) = master[key] + sum(w[key] for w in ws)
    S["memory_master"] = master
    S["memory_workers"] = ws
    S["memory_total"] = Dict(k => total(k) for k in ("rss_MB", "pss_MB", "peak_rss_MB"))
    @printf("memory: master RSS %.0f MB; workers RSS %s MB; total RSS %.0f MB, PSS %.0f MB, sum of peaks %.0f MB\n",
            master["rss_MB"], join([@sprintf("%.0f", w["rss_MB"]) for w in ws], "/"), total("rss_MB"),
            total("pss_MB"), total("peak_rss_MB"))
    mkpath(dirname(OUT))
    open(io -> JSON.json(io, S; pretty = true), OUT, "w")
    println("wrote $OUT")
    rmprocs(workers())
end

main()
