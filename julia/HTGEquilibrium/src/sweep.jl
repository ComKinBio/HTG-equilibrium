# Temperature sweeps and preflight checks (Python: htg/sweep.py, and Model.sweep
# and Model.check of htg/model.py).
#
# A sweep solves every point of the case's temperature-pressure grid, chain by
# chain (see `chains` and `solve_chain`): with continuation each point of a
# chain starts from the previous solution. With workers > 1 the chains are
# solved at the same time on threads, through a `WorkerPool` (queue.jl).

"""The progress line of one point, as Python prints it."""
point_text(T, P, outcome, start) =
    "$(py_format_g(T, 6)) °C, $(py_format_g(P, 6)) bar: " *
    "$(outcome isa Result ? "solved" : "FAILED: $outcome") (start: $start)"

"""
    run_sweep(model, case; progress=println, workers=1) -> SweepResult

Solve `case` at every point of its temperature-pressure grid.

With continuation each chain (rising temperature at one pressure, or the order
of explicit points) starts each point from the previous successful solution
when the selected (phase, species) pairs are unchanged; otherwise from a fresh
feasible LP. A failed point is recorded and the sweep continues. A heating path
can follow a local branch; smooth profiles do not prove global equilibrium.

`progress(text)` is called after every point (`nothing`: silent). workers > 1
solves independent chains at the same time on threads (start Julia with
several threads, e.g. `julia -t 4`); the model is shared.
"""
function run_sweep(model::Model, case::Case; progress=println, workers=1)
    progress = progress === nothing ? (_ -> nothing) : progress
    BLAS.set_num_threads(1)   # M1: one BLAS thread per solve (Ipopt/MUMPS are faster; no oversubscription)
    chain_list = chains(case)
    if workers > 1 && length(chain_list) > 1
        return _parallel(model, case, progress, min(workers, length(chain_list)))
    end
    outcomes = Dict{Tuple{Float64,Float64},Union{Result,String}}()
    on_point(T, P, outcome, start) = (outcomes[(T, P)] = outcome; progress(point_text(T, P, outcome, start)))
    for chain in chain_list
        solve_chain(model, case, chain, on_point)
    end
    return SweepResult(case, [(T, P, outcomes[(T, P)]) for (T, P) in grid(case)])
end

function _parallel(model, case, progress, workers)
    progress("Solving $(length(grid(case))) points with $workers workers (threads; " *
             "$(Threads.nthreads(:default)) available)")
    on_point(run, T, P, outcome, start) = progress(point_text(T, P, outcome, start))
    run = WorkerPool(workers; model, warm_up=false, on_point) do pool
        run = submit(pool, case)
        wait(run)
        run
    end
    return sweep_result(run)
end

"""
    sweep(model, case; progress=println, workers=1) -> SweepResult

Equilibrium at every point of the case (Python: `Model.sweep`); see `run_sweep`.
"""
sweep(model::Model, case::Case; progress=println, workers=1) = run_sweep(model, case; progress, workers)

"""
    preflight(model, case) -> DataTable

Check every point without the equation of state (columns T_C, P_bar, status
"ready" or "invalid", aqueous_removed, error).

Verifies species names, standard-state and SAFT data coverage at each
temperature, aqueous validity limits (with the case's water model), and
feasibility of the full set of balances and constraints. Passing does not
guarantee that the solver converges.
"""
function preflight(model::Model, case::Case)
    registry = model.registry
    rows = OrderedDict{String,Any}[]
    for (T, P_bar) in grid(case)
        point = at(case, T, P_bar)
        row = OrderedDict{String,Any}("T_C" => T, "P_bar" => P_bar, "status" => "ready", "aqueous_removed" => nothing,
                                      "error" => "")
        try
            phases, _, _, _, _, report = apply_phase_rules(point, registry)
            row["aqueous_removed"] = report["aqueous_disabled"]
            setup(model, point; build_phases=false)
            T_K, P = T + 273.15, point.pressure_bar * 1e5
            for (phase, names) in phases, name in names
                if phase == "solid"
                    solid_g0(model, [name], T_K, P)
                elseif !haskey(registry.ion, name)
                    nasa_g0(model, [name], T_K)
                    resolve_parameters(registry, model.databases, T_K, [name])
                end
            end
            if !isempty(phases["aqueous"]) && point.aqueous_model == "electrolyte"
                water = water_properties(water_model_of(point), T_K, P)
                if !(water.rho >= 350 && water.epsilon > 1)
                    rho = isnan(water.rho) ? "nan" : @sprintf("%.0f", water.rho)
                    throw(EquilibriumError("Aqueous electrolyte model invalid here: water density $rho kg/m3 " *
                                           "(needs >= 350). Set aqueous cutoff_C below this temperature."))
                end
            end
        catch exc
            exc isa InterruptException && rethrow()
            row["status"] = "invalid"
            row["error"] = "$(py_exception_name(exc)): $(_exception_text(exc))"
        end
        push!(rows, row)
    end
    return DataTable(rows)
end

"""`Model.check`: the preflight table of `case` (see `preflight`)."""
check(model::Model, case::Case) = preflight(model, case)

"""
    preflight_messages(model, case) -> Vector{String}

Readable reasons the case cannot run (empty: none found). First the
carbon-bound checks (contradictory or impossible bounds); if they pass, the
per-point preflight, as "T °C, P bar: error" for each failed point.
"""
function preflight_messages(model::Model, case::Case)
    solids = model.registry.solid
    carbon_solids = OrderedDict{String,Bool}(name => py_truth(get(solids[name]["elements"], "C", 0))
                                             for name in get(case.phases, "solid", String[]) if haskey(solids, name))
    messages = carbon_bound_messages(case, carbon_solids)
    isempty(messages) || return messages
    table = preflight(model, case)
    return ["$(py_format_g(table["T_C"][i], 6)) °C, $(py_format_g(table["P_bar"][i], 6)) bar: $(table["error"][i])"
            for i in 1:nrows(table) if table["status"][i] != "ready"]
end
