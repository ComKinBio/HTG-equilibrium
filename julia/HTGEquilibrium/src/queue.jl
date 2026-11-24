# Running cases on a pool of worker threads (Python: htg/runner.py, which uses
# worker processes).
#
# One Julia process (M1 decision: threads, not Distributed). Every submitted
# case is split into independent pieces (see `chains`): single points, or
# continuation chains that must run in order. All pieces of all runs go into one
# first-in, first-out queue; an idle worker (a task on a thread) takes the next
# piece. One large case spreads over all workers, and with several cases the one
# submitted first finishes first. The worker count limits how many chains are
# solved at the same time; the threads available are those Julia was started
# with (`julia -t N`, `JULIA_NUM_THREADS`).
#
#     pool = WorkerPool(4)                        # loads the model (data/) on start
#     run = submit(pool, case, "results/my_case")
#     wait(run)
#     summary(sweep_result(run))
#     shutdown(pool)
#
# Thread safety: the model is only read, so all workers share it; each point
# builds its own phase objects and Ipopt problem (solve). Ipopt serialises MUMPS
# with a mutex. BLAS runs on one thread (set on start). The pool's state (runs,
# queue, outcomes) is changed only under the pool's lock; the callbacks
# on_point and on_done are called under that lock, from a worker thread, one at
# a time (as Python calls them from its dispatcher thread).

const GB = 1024^3
"""Memory one more solver thread needs [GB] (M1 spike: ~0.1 GB per thread; with margin)."""
const WORKER_MEMORY_GB = 0.25

# ---------------------------------------------------------------- solving a chain
"""
    solve_chain(model, case, chain, on_point; stop=() -> false)

Solve the (T, P) points of `chain` in order; call `on_point(T, P, outcome,
start)` after each, where outcome is a `Result` or an error message (a
`SolveError`'s text). With continuation each point starts from the previous
successful solution when the selected (phase, species) pairs are unchanged;
otherwise from a fresh feasible LP. `stop()` is asked before each point (the
pool's shutdown); the remaining points are then left out.
"""
function solve_chain(model::Model, case::Case, chain, on_point; stop=() -> false)
    previous = nothing
    for (T, P) in chain
        stop() && return nothing
        point = at(case, T, P)
        start = "feasible LP"
        if py_truth(case.continuation) && previous !== nothing
            phases = apply_phase_rules(point, model.registry)[1]
            signature = Set(Pair2[(p, s) for (p, names) in phases for s in names])
            if Set(keys(previous.amounts)) == signature
                point.initial_amounts = OrderedDict{Pair2,Float64}(previous.amounts)
                start = "solution at $(py_format_g(previous.temperature_C, 6)) °C, " *
                        "$(py_format_g(previous.pressure_bar, 6)) bar"
            else
                point.initial_amounts = nothing
                start = "feasible LP (phase selection changed)"
            end
        end
        outcome = try
            previous = solve(model, point)
        catch exc
            exc isa SolveError || rethrow()
            exc.msg
        end
        on_point(T, P, outcome, start)
    end
    return nothing
end

"""
    default_workers() -> Int

Number of workers: min(threads Julia was started with, CPU threads - 1,
available memory / WORKER_MEMORY_GB), at least 1.
"""
function default_workers()
    cores = max(1, Sys.CPU_THREADS - 1)
    memory = floor(Int, Sys.free_memory() / (WORKER_MEMORY_GB * GB))
    return max(1, min(Threads.nthreads(:default), cores, memory))
end

# ---------------------------------------------------------------- runs
"""
    Run

A submitted case; updated by the pool while it runs. Fields as in Python:
`id`, `case`, `output_dir`, `total` (points), `status` ("queued", "running",
"cancelling", "done", "cancelled"), `outcomes` ((T, P) => (Result or error
message, start)), `log`, `submitted`, `started`, `finished` (times),
`tasks_left`; and `model`, the model the run uses. Properties as in Python:
`run.name`, `run.done_points`, `run.failed_points`, `run.finished_running`.
`wait(run, timeout)`, `sweep_result(run)`.
"""
mutable struct Run
    id::Int
    case::Case
    output_dir::Union{Nothing,String}
    total::Int
    status::String
    outcomes::OrderedDict{Tuple{Float64,Float64},Tuple{Union{Result,String},String}}
    log::Vector{String}
    submitted::Float64
    started::Union{Nothing,Float64}
    finished::Union{Nothing,Float64}
    tasks_left::Int
    model::Model
    lock::ReentrantLock          # the pool's lock
    done_event::Base.Event
end

function Base.getproperty(run::Run, name::Symbol)
    name === :name && return getfield(run, :case).name
    name === :done_points && return lock(() -> length(getfield(run, :outcomes)), getfield(run, :lock))
    name === :failed_points &&
        return lock(() -> count(o -> o[1] isa String, values(getfield(run, :outcomes))), getfield(run, :lock))
    name === :finished_running && return getfield(run, :status) in ("done", "cancelled")
    return getfield(run, name)
end
Base.propertynames(::Run) = (fieldnames(Run)..., :name, :done_points, :failed_points, :finished_running)
Base.show(io::IO, run::Run) = print(io, "Run(#$(run.id) $(run.name): $(run.status), $(run.done_points)/$(run.total) points)")

"""
    wait(run, timeout=nothing) -> Bool

Wait until the run has finished (done or cancelled); false if `timeout` seconds passed first.
"""
function Base.wait(run::Run, timeout=nothing)
    timeout === nothing && (wait(run.done_event); return true)
    return timedwait(() -> run.finished_running, Float64(timeout); pollint=0.01) === :ok
end

"""Results so far, in grid order (points not solved yet are left out)."""
function sweep_result(run::Run)
    points = lock(run.lock) do
        [(T, P, run.outcomes[(T, P)][1]) for (T, P) in grid(run.case) if haskey(run.outcomes, (T, P))]
    end
    return SweepResult(run.case, points)
end

"""
    live_tables(run) -> RunTables

Tables of a run in the pool (finished or still running), with the choices saved
in its output folder (Python: `htg.gui.runs.live_tables`).
"""
function live_tables(run::Run)
    sweep = sweep_result(run)
    comparison = nothing
    if py_truth(run.case.comparison) && isfile(String(run.case.comparison["file"]))
        comparison = read_csv(String(run.case.comparison["file"]))
    end
    folder = run.output_dir
    choices = folder === nothing ? Dict() : read_choices(folder)
    return RunTables("$(run.name) (this session, #$(run.id))", run.case, summary(sweep), species_table(sweep),
                     phase_summary(sweep); comparison, path=folder, alternatives=alternatives_table(sweep), choices)
end

"""A piece of a run: one chain of points, solved in order by one worker."""
mutable struct ChainTask
    id::Int
    run::Run
    chain::Vector{Tuple{Float64,Float64}}
end

"""The first point of a task's chain without an outcome yet (the one being solved)."""
_current_point(::Nothing) = nothing
function _current_point(task::ChainTask)
    i = findfirst(p -> !haskey(task.run.outcomes, p), task.chain)
    return i === nothing ? nothing : task.chain[i]
end

mutable struct PoolWorker
    id::Int
    state::String   # starting, warming, idle, busy, stopped, crashed
    task::Union{Nothing,ChainTask}
    retire::Bool
    message::String
    warm_up::Bool
    handle::Union{Nothing,Task}
end

# ---------------------------------------------------------------- pool
"""
    WorkerPool(workers=default_workers(); model=nothing, model_args=nothing, warm_up=true,
               on_point=nothing, on_done=nothing, chain_solver=solve_chain)

Worker tasks on threads with one FIFO queue of chains (Python: `WorkerPool`,
with processes).

- `workers`: number of chains solved at the same time
- `model`: the `Model` to use; else it is loaded from `model_args` (the
  arguments of `load_model` as a dict: "data_dir", "species_file",
  "nasa_names", "reference_temperature_C", "user_data"; Python's
  `Model.load_args`) when the pool starts
- `warm_up`: one worker solves a small case first so that the code is compiled
  before the first real point (all threads share the compiled code)
- `on_point(run, T, P, outcome, start)`, `on_done(run)`: optional callbacks,
  called from a worker thread under the pool's lock
- `chain_solver(model, case, chain, on_point; stop)`: solves one chain
  (tests replace it with a stand-in)

`WorkerPool(f, workers; ...)` (do-block) starts the pool, calls `f(pool)` and
shuts it down (Python: `with WorkerPool(...) as pool`).
"""
mutable struct WorkerPool
    target::Int
    model::Union{Nothing,Model}
    model_args::OrderedDict{String,Any}
    warm_up::Bool
    on_point::Any
    on_done::Any
    chain_solver::Any
    lock::ReentrantLock
    cond::Threads.Condition
    workers::OrderedDict{Int,PoolWorker}
    pending::Vector{ChainTask}
    runs::Vector{Run}
    last_id::Int
    started::Bool
    stopping::Bool
    warmed::Bool
end

function WorkerPool(workers=nothing; model=nothing, model_args=nothing, warm_up=true, on_point=nothing,
                    on_done=nothing, chain_solver=solve_chain)
    lk = ReentrantLock()
    args = OrderedDict{String,Any}(String(k) => v for (k, v) in something(model_args, Dict()))
    if model !== nothing && isempty(args)
        args = copy(model.load_args)
    end
    return WorkerPool(Int(something(workers, default_workers())), model, args, warm_up, on_point, on_done,
                      chain_solver, lk, Threads.Condition(lk), OrderedDict{Int,PoolWorker}(), ChainTask[], Run[],
                      0, false, false, false)
end

function WorkerPool(f::Function, workers=nothing; kwargs...)
    pool = start(WorkerPool(workers; kwargs...))
    try
        return f(pool)
    finally
        shutdown(pool)
    end
end

_next_id!(pool::WorkerPool) = (pool.last_id += 1)

_arg(args, key, default=nothing) = (v = get(args, key, nothing); v === nothing ? default : v)

"""Load a model from Python-style `Model.load_args`."""
load_model_from_args(args) =
    load_model(_arg(args, "data_dir"); species_file=_arg(args, "species_file"), nasa_names=_arg(args, "nasa_names"),
               reference_temperature_C=_arg(args, "reference_temperature_C", 350.0), user_data=_arg(args, "user_data"))

function _ensure_model!(pool::WorkerPool)
    pool.model === nothing || return pool.model
    model = load_model_from_args(pool.model_args)
    lock(pool.lock) do
        pool.model === nothing && (pool.model = model)
    end
    return pool.model
end

"""Start the workers (and load the model, if not given). Returns the pool."""
function start(pool::WorkerPool)
    pool.started && return pool
    _ensure_model!(pool)
    BLAS.set_num_threads(1)   # M1: one BLAS thread per solve
    lock(pool.lock) do
        if !pool.started
            pool.started = true
            for _ in 1:pool.target
                _spawn!(pool)
            end
        end
    end
    return pool
end

"""
    submit(pool, case, output_dir=nothing) -> Run

Queue `case`. Results are saved to `output_dir` when it finishes (summary.csv,
species.csv, phases.csv, alternatives.csv, case.toml, comparison.csv). The run
uses the pool's current model.
"""
function submit(pool::WorkerPool, case::Case, output_dir=nothing)
    validate(case)
    chain_list = chains(case)
    model = _ensure_model!(pool)
    run = lock(pool.lock) do
        run = Run(_next_id!(pool), case, output_dir === nothing ? nothing : String(output_dir),
                  sum(length, chain_list; init=0), "queued",
                  OrderedDict{Tuple{Float64,Float64},Tuple{Union{Result,String},String}}(), String[], time(),
                  nothing, nothing, length(chain_list), model, pool.lock, Base.Event())
        for chain in chain_list
            push!(pool.pending, ChainTask(_next_id!(pool), run, chain))
        end
        push!(pool.runs, run)
        isempty(chain_list) && _finish!(pool, run, "done")
        notify(pool.cond; all=true)
        run
    end
    start(pool)
    return run
end

"""Drop the run's queued pieces; pieces already running finish their current chain."""
function cancel(pool::WorkerPool, run::Run)
    lock(pool.lock) do
        run.finished_running && return nothing
        kept = ChainTask[]
        for task in pool.pending
            if task.run === run
                run.tasks_left -= 1
            else
                push!(kept, task)
            end
        end
        pool.pending = kept
        push!(run.log, "cancelled")
        run.status = "cancelling"
        run.tasks_left <= 0 && _finish!(pool, run, "cancelled")
        return nothing
    end
end

_active(w::PoolWorker) = !(w.state in ("stopped", "crashed")) && !w.retire

"""Change the number of workers (idle extra workers stop; busy ones after their chain)."""
function set_workers(pool::WorkerPool, n)
    lock(pool.lock) do
        pool.target = max(1, Int(n))
        active = [w for w in values(pool.workers) if _active(w)]
        if pool.started
            for _ in 1:(pool.target - length(active))
                _spawn!(pool)
            end
        end
        for w in active[pool.target+1:end]
            w.retire = true
        end
        notify(pool.cond; all=true)
    end
    return pool
end

"""
    restart(pool, model_args) -> Bool
    restart(pool, model::Model) -> Bool

Use a new model (after the user's species data changed) for runs submitted from
now on (Python: `WorkerPool.restart`, which replaced the worker processes).
Running work is never stopped and a run never mixes old and new data: runs
already submitted finish with the model they started with. Returns true when no
run was unfinished (everything uses the new model now), false otherwise.
"""
restart(pool::WorkerPool, model_args) = restart(pool, load_model_from_args(
    OrderedDict{String,Any}(String(k) => v for (k, v) in something(model_args, Dict()))))
function restart(pool::WorkerPool, model::Model)
    lock(pool.lock) do
        pool.model = model
        pool.model_args = copy(model.load_args)
        return !any(run -> !run.finished_running, pool.runs)
    end
end

"""
    status(pool) -> OrderedDict

Worker states and queue length (for display): "workers" (each with "id",
"state", "message", "case" and "point", the (T_C, P_bar) it is solving or
nothing), "queued_pieces", "target".
"""
function status(pool::WorkerPool)
    lock(pool.lock) do
        workers = [OrderedDict{String,Any}("id" => w.id, "state" => w.state, "message" => w.message,
                                           "case" => w.task === nothing ? "" : w.task.run.name,
                                           "point" => _current_point(w.task))
                   for w in values(pool.workers) if w.state != "stopped"]
        return OrderedDict{String,Any}("workers" => workers, "queued_pieces" => length(pool.pending),
                                       "target" => pool.target)
    end
end

"""
    shutdown(pool; wait=true)

Stop the workers: each finishes the point it is solving; the rest of its chain
and the queued pieces are not solved (their runs stay unfinished, as in Python).
With `wait`, return when every worker has stopped.
"""
function shutdown(pool::WorkerPool; wait=true)
    handles = lock(pool.lock) do
        pool.stopping = true
        notify(pool.cond; all=true)
        [w.handle for w in values(pool.workers) if w.handle !== nothing]
    end
    if wait
        for handle in handles
            Base.wait(handle)
        end
    end
    return nothing
end

# ---------------------------------------------------------------- internals (called under the lock)
function _spawn!(pool::WorkerPool)
    worker = PoolWorker(_next_id!(pool), "starting", nothing, false, "", pool.warm_up && !pool.warmed, nothing)
    pool.warmed = pool.warmed || worker.warm_up
    pool.workers[worker.id] = worker
    worker.handle = Threads.@spawn :default _worker_loop(pool, worker)
    return worker
end

"""
The commented case-file template (`htg new`; the warm-up solve): a copy of the
Python app's `src/htg/templates/case_template.toml` (kept identical, tested).
"""
const CASE_TEMPLATE = normpath(joinpath(@__DIR__, "..", "templates", "case_template.toml"))

function _warm_up(pool::WorkerPool, worker::PoolWorker)
    lock(pool.lock) do
        worker.state, worker.message = "warming", "compiling the equation of state"
    end
    try
        solve(pool.model, at(load_case(CASE_TEMPLATE), 300.0, 240.0))
    catch exc
        exc isa InterruptException && rethrow()
    end
end

function _worker_loop(pool::WorkerPool, worker::PoolWorker)
    try
        worker.warm_up && _warm_up(pool, worker)
        while true
            task = lock(pool.lock) do
                _next_task!(pool, worker)
            end
            task === nothing && return nothing
            _run_task(pool, worker, task)
        end
    catch exc
        text = sprint(showerror, exc, catch_backtrace())
        lock(pool.lock) do
            worker.state, worker.message = "crashed", first(text, 500)
            task, worker.task = worker.task, nothing
            if task !== nothing
                _fail_rest!(task, "worker stopped unexpectedly")
                _task_finished!(pool, task)
            end
            # replace the worker, as Python replaces a crashed process
            pool.stopping || worker.retire || _spawn!(pool)
        end
    end
    return nothing
end

"""The next piece of the queue for `worker` (waits while the queue is empty); nothing: stop."""
function _next_task!(pool::WorkerPool, worker::PoolWorker)
    while true
        if pool.stopping || worker.retire
            worker.state, worker.task = "stopped", nothing
            return nothing
        end
        if !isempty(pool.pending)
            task = popfirst!(pool.pending)
            run = task.run
            if run.status == "queued"
                run.status, run.started = "running", time()
            end
            worker.task, worker.state, worker.message = task, "busy", ""
            return task
        end
        worker.state = "idle"
        wait(pool.cond)
    end
end

function _run_task(pool::WorkerPool, worker::PoolWorker, task::ChainTask)
    run = task.run
    on_point(T, P, outcome, start) = lock(() -> _record_point!(pool, task, T, P, outcome, start), pool.lock)
    try
        pool.chain_solver(run.model, run.case, task.chain, on_point; stop=() -> pool.stopping)
    catch exc
        exc isa InterruptException && rethrow()
        text = sprint(showerror, exc, catch_backtrace())
        lock(() -> _fail_rest!(task, "internal error: $(first(text, 800))"), pool.lock)
    end
    lock(pool.lock) do
        worker.task = nothing
        worker.state = "idle"
        complete = all(p -> haskey(run.outcomes, p), task.chain)
        (complete || !pool.stopping) && _task_finished!(pool, task)
    end
    return nothing
end

function _record_point!(pool::WorkerPool, task::ChainTask, T, P, outcome, start)
    run = task.run
    run.outcomes[(Float64(T), Float64(P))] = (outcome, String(start))
    push!(run.log, point_text(T, P, outcome, start))
    pool.on_point === nothing || _callback(pool.on_point, run, T, P, outcome, start)
    return nothing
end

function _fail_rest!(task::ChainTask, reason)
    for (T, P) in task.chain
        if !haskey(task.run.outcomes, (T, P))
            task.run.outcomes[(T, P)] = (reason, "")
            push!(task.run.log, "$(py_format_g(T, 6)) °C, $(py_format_g(P, 6)) bar: FAILED: $reason")
        end
    end
end

function _task_finished!(pool::WorkerPool, task::ChainTask)
    run = task.run
    run.tasks_left -= 1
    if run.tasks_left <= 0 && !run.finished_running
        _finish!(pool, run, run.status == "cancelling" ? "cancelled" : "done")
    end
end

function _finish!(pool::WorkerPool, run::Run, status)
    run.status, run.finished = status, time()
    if run.output_dir !== nothing && !isempty(run.outcomes)
        try
            save(sweep_result(run), run.output_dir)
            push!(run.log, "results saved to $(run.output_dir)")
        catch exc
            exc isa InterruptException && rethrow()
            push!(run.log, "saving results failed: $(sprint(showerror, exc))")
        end
    end
    pool.on_done === nothing || _callback(pool.on_done, run)
    notify(run.done_event)
end

function _callback(f, args...)
    try
        f(args...)
    catch exc
        exc isa InterruptException && rethrow()
        @error "callback failed" exception = (exc, catch_backtrace())
    end
end
