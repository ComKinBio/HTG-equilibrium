# What the tabs share (Python: htg/gui/app.py State, gui/common.py
# worker_summary, gui/runs.py output_folder).

"""
    worker_summary(status) -> String

The status line text for `status(pool)` (Python: `gui/common.py worker_summary`):
"Workers: 3 busy, 1 idle · 2 piece(s) waiting".
"""
function worker_summary(status)
    counts = Dict{String,Int}()
    for worker in status["workers"]
        counts[worker["state"]] = get(counts, worker["state"], 0) + 1
    end
    order = ("busy", "idle", "warming", "starting", "stopping", "crashed")
    parts = ["$(counts[s]) $s" for s in order if get(counts, s, 0) > 0]
    text = "Workers: " * (isempty(parts) ? "none" : join(parts, ", "))
    queued = status["queued_pieces"]
    queued > 0 && (text *= " · $queued piece(s) waiting")
    return text
end

"""results/<case name>-<date-time>, so runs never overwrite each other (Python: `gui/runs.py output_folder`)."""
output_folder(results_dir, case_name, now=Dates.now()) =
    joinpath(results_dir, "$(case_name)-$(Dates.format(now, "yyyymmdd-HHMMSS"))")

"""Python's `Path(p).resolve()`: absolute, with the links of the existing part resolved."""
function resolve_path(p)
    path = abspath(expanduser(String(p)))
    existing, rest = path, String[]
    while !ispath(existing) && dirname(existing) != existing
        pushfirst!(rest, basename(existing))
        existing = dirname(existing)
    end
    ispath(existing) && (existing = realpath(existing))
    return isempty(rest) ? existing : joinpath(existing, rest...)
end

"""
    State(model=nothing, pool=nothing; results_dir="results", cases_dir="cases", settings=nothing)

Everything shared by the tabs (one local user; Python: `gui/app.py State`).

- `model`, `pool`: the model and the worker pool. The app starts with neither
  and fills them in from a background task (`load!`), so the window opens
  first; `ready(state)` says when both are there, `problem` why loading failed.
- `settings`: app preferences (`HTGEquilibrium.Settings`), read once here,
  written by the tabs with `settings[key] = value`. Default: `load_settings()`.
- `results_dir`, `cases_dir`: absolute paths.
- `phase_options` ("gas", "aqueous", "organic_liquid", "solid" => species
  names) and `feed_species`: follow the model (`use_model!`).
- `document`: the case being edited (`CaseDocument`; a new case at start).

Not ported yet (they come with the tabs that use them): `sources`, `run`, `tables`.
"""
mutable struct State
    model::Union{Nothing,Model}
    pool::Union{Nothing,WorkerPool}
    settings::Settings
    results_dir::String
    cases_dir::String
    phase_options::OrderedDict{String,Vector{String}}
    feed_species::Vector{String}
    problem::String
    document::CaseDocument
end

function State(model=nothing, pool=nothing; results_dir="results", cases_dir="cases", settings=nothing)
    state = State(nothing, pool, settings === nothing ? load_settings() : settings, resolve_path(results_dir),
                  resolve_path(cases_dir), OrderedDict{String,Vector{String}}(), String[], "", new_document())
    model === nothing || use_model!(state, model)
    return state
end

"""Use `model` (for example after the user species changed): the species lists follow it."""
function use_model!(state::State, model::Model)
    catalogue = species(model)
    names(phase) = [r.Species for r in catalogue if r.Phase == phase]
    molecular, ions, solids = names("molecular"), names("aqueous ion"), names("solid")
    state.phase_options = OrderedDict("gas" => molecular, "aqueous" => vcat(molecular, ions),
                                      "organic_liquid" => molecular, "solid" => solids)
    state.feed_species = sort!(collect(union(keys(all_molecular(model.registry)), keys(model.registry.solid))))
    state.model = model
    return state
end

"""True when the model is loaded and the worker pool is running."""
ready(state::State) = state.model !== nothing && state.pool !== nothing

"""CSV files below the cases folder: (path relative to `base_dir`, path shown) (Python: State.comparison_files)."""
function comparison_files(state::State, base_dir)
    isdir(state.cases_dir) || return Tuple{String,String}[]
    files = String[]
    for (dir, _, names) in walkdir(state.cases_dir), name in names
        endswith(name, ".csv") && push!(files, joinpath(dir, name))
    end
    return [(replace(relpath(f, base_dir), '\\' => '/'), relpath(f, state.cases_dir)) for f in sort!(files)]
end

"""Queue `case`; its results go to a new folder below `results_dir`."""
function HTGEquilibrium.submit(state::State, case::Case)
    ready(state) || throw(ArgumentError("The app is still loading the species data. Try again in a moment."))
    return submit(state.pool, case, output_folder(state.results_dir, case.name))
end

"""The text of the status line in the header."""
function status_text(state::State)
    isempty(state.problem) || return "Could not load the species data: $(state.problem)"
    state.pool === nothing && return "Loading the species data"
    return worker_summary(status(state.pool))
end

"""
    load!(state; workers=nothing, model_args=Dict(), warm_up=true, out=stdout) -> state

Load the model and start the worker pool (Python did this before opening the
window; the app now runs it in a background task). `workers`: worker threads
(default `default_workers()`; at most the threads Julia has). Errors are kept
in `state.problem` and printed.
"""
function load!(state::State; workers=nothing, model_args=Dict{String,Any}(), warm_up=true, out::IO=stdout)
    try
        model = load_model_from_args(OrderedDict{String,Any}(String(k) => v for (k, v) in model_args))
        use_model!(state, model)
        n = workers === nothing ? default_workers() : worker_threads(workers, out)
        pool = WorkerPool(n; model, warm_up)
        println(out, "Starting $n worker thread(s) (the first results take about a minute while Julia compiles)")
        flush(out)
        HTGEquilibrium.start(pool)
        state.pool = pool
    catch exc
        exc isa InterruptException && rethrow()
        state.problem = _one_line(exc)
        println(stderr, "Error: could not load the species data: ", sprint(showerror, exc))
        flush(stderr)
    end
    return state
end

_one_line(exc) = first(join(split(strip(sprint(showerror, exc)), '\n'; keepempty=false), " "), 300)
