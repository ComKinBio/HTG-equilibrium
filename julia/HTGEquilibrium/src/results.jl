# Results of equilibrium calculations, as objects and tables (Python: htg/results.py).
#
# Amounts are in mol per mol of dry-feed carbon. Tables also give yields per kg
# of original dry feed (salts are extra mass and not included in the basis).
# Tables are `DataTable`s (tables.jl) and are written as pandas writes them, so
# saved runs from either version open in either version.

"""(phase, species) pair, the key of amounts and species data."""
const Pair2 = Tuple{String,String}

"""
    Alternative(label, gibbs, amounts, n_starts, start_labels; phase_gibbs=OrderedDict())

One distinct result at a point (the converged starts that found it).

    label         "A", "B", ... by increasing Gibbs energy; A is the main result
    gibbs         total Gibbs energy [J/mol dry-feed C]
    amounts       {(phase, species) => mol / mol dry-feed C}
    n_starts      number of starts that found it; start_labels: their names
    phase_gibbs   {phase => Gibbs energy [J/mol dry-feed C]}
"""
mutable struct Alternative
    label::String
    gibbs::Float64
    amounts::OrderedDict{Pair2,Float64}
    n_starts::Int
    start_labels::Vector{String}
    phase_gibbs::OrderedDict{String,Float64}
end
Alternative(label, gibbs, amounts, n_starts, start_labels; phase_gibbs=OrderedDict{String,Float64}()) =
    Alternative(String(label), Float64(gibbs), OrderedDict{Pair2,Float64}(amounts), Int(n_starts),
                collect(String, start_labels), OrderedDict{String,Float64}(phase_gibbs))
Alternative(label, gibbs, amounts, n_starts, start_labels, phase_gibbs::AbstractDict) =
    Alternative(label, gibbs, amounts, n_starts, start_labels; phase_gibbs)

"""
    Result

Equilibrium at one temperature and pressure (Python: `htg.results.Result`, same fields).

    amounts                {(phase, species) => mol / mol dry-feed C}
    gibbs_J_per_mol_C      total Gibbs energy at the solution
    gas_carbon_fraction    carbon in gas / all feed carbon (including salt carbon)
    carbon_bound_violation largest violation of the carbon bounds [fraction of feed carbon]
    phase_rules            what the aqueous cut-off did at this temperature
                           ("dropped_constraints": bounds dropped at the cut-off)
    optimizer              solver diagnostics (KKT residuals, starts, results, ...)
    alternatives           the distinct results [Alternative], A first (A = this result)
    water_model            the case's water model name (Julia only; "SUPCRT92" or "IAPWS-95")
"""
mutable struct Result
    temperature_C::Float64
    pressure_bar::Float64
    amounts::OrderedDict{Pair2,Float64}
    initial_amounts::OrderedDict{Pair2,Float64}
    gibbs_J_per_mol_C::Float64
    phase_gibbs_J_per_mol_C::OrderedDict{String,Float64}
    balance_residual::Float64
    charge_residual::Float64
    iterations::Int
    eos_trial_errors::Vector{String}
    feed_C_mol_per_g_dry::Float64
    elements::Vector{String}
    feed_atoms::Vector{Float64}
    product_atoms::Vector{Float64}
    gas_carbon_fraction::Float64
    carbon_bound_violation::Float64
    phase_rules::OrderedDict{String,Any}
    optimizer::OrderedDict{String,Any}
    atomic_weights::Dict{String,Float64}
    species_data::OrderedDict{Pair2,OrderedDict{String,Any}}
    alternatives::Vector{Alternative}
    water_model::String
end

const RESULT_FIELDS = fieldnames(Result)

"""Result(; fields...) with Python's defaults (alternatives = [], water_model = "SUPCRT92")."""
function Result(; alternatives=Alternative[], water_model="SUPCRT92", kwargs...)
    for k in keys(kwargs)
        k in RESULT_FIELDS || throw(ArgumentError("Result.__init__() got an unexpected keyword argument $(py_repr(string(k)))"))
    end
    values = Any[]
    for f in RESULT_FIELDS
        f === :alternatives ? push!(values, collect(Alternative, alternatives)) :
        f === :water_model ? push!(values, String(water_model)) :
        haskey(kwargs, f) ? push!(values, kwargs[f]) :
        throw(ArgumentError("Result.__init__() missing 1 required positional argument: $(py_repr(string(f)))"))
    end
    return Result(values...)
end

"""Number of distinct results at the point (at least 1)."""
n_results(r::Result) = max(1, length(r.alternatives))

"""Short notes for the summary table, e.g. "aqueous carbon bound dropped"."""
notes(r::Result) = join(("$name bound dropped" for name in get(r.phase_rules, "dropped_constraints", String[])), "; ")

"""
    for_alternative(result, label) -> Result

This result with the amounts and Gibbs energies of alternative `label`.
"""
function for_alternative(r::Result, label)
    label == "A" && !any(a -> a.label == "A", r.alternatives) && return r
    k = findfirst(a -> a.label == label, r.alternatives)
    k === nothing && throw(KeyError("No result $(py_repr(String(label))) at $(py_format_g(r.temperature_C, 6)) °C, " *
                                    "$(py_format_g(r.pressure_bar, 6)) bar."))
    alt = r.alternatives[k]
    amounts = OrderedDict{Pair2,Float64}(pair => Float64(get(alt.amounts, pair, 0.0)) for pair in keys(r.amounts))
    feed_C = r.feed_atoms[findfirst(==("C"), r.elements)]
    gas_C = sum((v * get(r.species_data[pair]["elements"], "C", 0) for (pair, v) in amounts if pair[1] == "gas");
                init=0.0)
    product = Float64[sum((v * get(r.species_data[pair]["elements"], e, 0) for (pair, v) in amounts); init=0.0)
                      for e in r.elements]
    phase_gibbs = OrderedDict{String,Float64}(phase => 0.0 for phase in keys(r.phase_gibbs_J_per_mol_C))
    for (k, v) in alt.phase_gibbs
        phase_gibbs[k] = Float64(v)
    end
    out = deepcopy(r)
    out.amounts = amounts
    out.gibbs_J_per_mol_C = alt.gibbs
    out.phase_gibbs_J_per_mol_C = phase_gibbs
    out.product_atoms = product
    out.gas_carbon_fraction = gas_C / feed_C
    return out
end

# Python's __repr__: the point, G, gas carbon and the 12 largest amounts.
function Base.show(io::IO, ::MIME"text/plain", r::Result)
    major = sort([(v, p, s) for ((p, s), v) in r.amounts if v > 1e-3]; rev=true)
    print(io, "Result at $(py_format_g(r.temperature_C, 6)) °C, $(py_format_g(r.pressure_bar, 6)) bar: ",
          "G = $(py_format_g(r.gibbs_J_per_mol_C, 6)) J/mol C, ",
          "$(@sprintf("%.1f", 100 * r.gas_carbon_fraction)) % of carbon in gas")
    for (v, p, s) in first(major, 12)
        print(io, "\n  ", rpad(p, 15), " ", rpad(s, 20), " ", py_format_g(v, 4), " mol/mol C")
    end
end
Base.show(io::IO, r::Result) =
    print(io, "Result($(py_format_g(r.temperature_C, 6)) °C, $(py_format_g(r.pressure_bar, 6)) bar, ",
          "G = $(py_format_g(r.gibbs_J_per_mol_C, 6)) J/mol C)")

# ---------------------------------------------------------------- tables of one result
"""File of a saved run that holds the chosen results."""
const CHOICES_FILE = "choices.toml"
"""Columns of alternatives.csv that are not species-table columns."""
const ALTERNATIVE_COLUMNS = ("result", "gibbs_J_per_mol_C", "n_starts", "phase_gibbs_J_per_mol_C", "start_labels")

_mw(aw, elements) = sum((aw[e] * c for (e, c) in elements); init=0.0)

"""
    species_table(result) -> DataTable

One row per (phase, species): amounts, yields per kg dry feed, carbon, molality
(columns T_C, P_bar, phase, group, species, formula, mol_per_mol_C,
mol_per_kg_dry_feed, g_per_kg_dry_feed, carbon_mol_per_kg_dry_feed,
molality_mol_per_kg_water).
"""
function species_table(r::Result)
    factor = r.feed_C_mol_per_g_dry * 1000.0  # mol/mol C -> mol/kg dry feed
    aw = r.atomic_weights
    water_kg = get(r.amounts, ("aqueous", "water"), 0.0) * (2 * aw["H"] + aw["O"]) / 1000.0
    rows = OrderedDict{String,Any}[]
    for ((phase, species), amount) in r.amounts
        data = r.species_data[(phase, species)]
        is_ion = py_truth(data["charge"])
        mw = _mw(aw, data["elements"])
        push!(rows, OrderedDict{String,Any}(
            "T_C" => r.temperature_C, "P_bar" => r.pressure_bar, "phase" => phase,
            "group" => is_ion ? "aqueous ions" : phase == "aqueous" ? "aqueous neutrals" : phase,
            "species" => species, "formula" => data["formula"],
            "mol_per_mol_C" => amount, "mol_per_kg_dry_feed" => amount * factor,
            "g_per_kg_dry_feed" => amount * factor * mw,
            "carbon_mol_per_kg_dry_feed" => amount * factor * Float64(get(data["elements"], "C", 0.0)),
            "molality_mol_per_kg_water" => is_ion && water_kg > 0 ? amount / water_kg : NaN))
    end
    return DataTable(rows)
end


"""
    phase_summary(result) -> DataTable

Per phase (all four, in `PHASES` order): total amount, mass, carbon share of all
feed carbon [%], Gibbs energy.
"""
function phase_summary(r::Result)
    feed_C = r.feed_atoms[findfirst(==("C"), r.elements)]
    factor = r.feed_C_mol_per_g_dry * 1000.0
    aw = r.atomic_weights
    rows = OrderedDict{String,Any}[]
    for phase in PHASES
        part = [(pair, v) for (pair, v) in r.amounts if pair[1] == phase]
        carbon = sum((v * get(r.species_data[pair]["elements"], "C", 0) for (pair, v) in part); init=0.0)
        mass = sum((v * factor * _mw(aw, r.species_data[pair]["elements"]) for (pair, v) in part); init=0.0)
        push!(rows, OrderedDict{String,Any}(
            "T_C" => r.temperature_C, "P_bar" => r.pressure_bar, "phase" => phase,
            "mol_per_mol_C" => sum((v for (_, v) in part); init=0.0), "g_per_kg_dry_feed" => mass,
            "carbon_percent" => 100 * carbon / feed_C,
            "gibbs_J_per_mol_C" => get(r.phase_gibbs_J_per_mol_C, phase, 0.0)))
    end
    return DataTable(rows)
end

"""
    summary(result) -> OrderedDict

The summary-table row of a solved point (Python: `Result.summary()`; extends `Base.summary`).
"""
Base.summary(r::Result) = OrderedDict{String,Any}(
    "T_C" => r.temperature_C, "P_bar" => r.pressure_bar, "status" => "solved",
    "gibbs_J_per_mol_C" => r.gibbs_J_per_mol_C,
    "gas_carbon_percent" => 100 * r.gas_carbon_fraction,
    "balance_residual" => r.balance_residual, "charge_residual" => r.charge_residual,
    "iterations" => r.iterations,
    "aqueous_removed" => get(r.phase_rules, "aqueous_disabled", nothing),
    "solver" => get(r.optimizer, "method", nothing),
    "descent_RT" => get(r.optimizer, "descent_RT", NaN),
    "n_results" => n_results(r), "carbon_bound_violation" => r.carbon_bound_violation,
    "notes" => notes(r), "error" => "")

# ---------------------------------------------------------------- sweeps
"""
    SweepResult(case, points; user_species=nothing)

Results over a temperature-pressure grid. `points`: [(T_C, P_bar, Result or
error message)] in grid order. Tables: `summary`, `species_table`,
`phase_summary`, `alternatives_table`; `HTGEquilibrium.save(sweep, directory)`
writes them.
"""
mutable struct SweepResult
    case::Case
    points::Vector{Tuple{Float64,Float64,Union{Result,String}}}
    user_species::Any   # Python: the user species run record (not written by Julia yet)
end
SweepResult(case::Case, points; user_species=nothing) =
    SweepResult(case, Tuple{Float64,Float64,Union{Result,String}}[(Float64(T), Float64(P), r) for (T, P, r) in points],
                user_species)

"""The solved points' results, in grid order."""
results(s::SweepResult) = Result[r for (_, _, r) in s.points if r isa Result]

"""One row per point: the result's summary, or status "failed" with the error."""
Base.summary(s::SweepResult) = DataTable(OrderedDict{String,Any}[r isa Result ? summary(r) :
                                          OrderedDict{String,Any}("T_C" => T, "P_bar" => P, "status" => "failed",
                                                                  "error" => r)
                                          for (T, P, r) in s.points])

function species_table(s::SweepResult)
    tables = [species_table(r) for r in results(s)]
    return isempty(tables) ? DataTable() : concat(tables)
end

function phase_summary(s::SweepResult)
    tables = [phase_summary(r) for r in results(s)]
    return isempty(tables) ? DataTable() : concat(tables)
end

const ALTERNATIVES_FIRST = ["T_C", "P_bar", "result", "gibbs_J_per_mol_C", "n_starts", "phase", "species",
                            "mol_per_mol_C"]

"""
    alternatives_table(sweep) -> DataTable

Long table of the distinct results at points with more than one result.
Columns: T_C, P_bar, result, gibbs_J_per_mol_C, n_starts, phase, species,
mol_per_mol_C, then the other species-table columns, the Gibbs energy of the
row's phase and the starts that found the result.
"""
function alternatives_table(s::SweepResult)
    tables = DataTable[]
    for r in results(s)
        length(r.alternatives) < 2 && continue
        for alt in r.alternatives
            table = species_table(for_alternative(r, alt.label))
            insert_column!(table, 2, "result", alt.label)
            insert_column!(table, 3, "gibbs_J_per_mol_C", Float64(alt.gibbs))
            insert_column!(table, 4, "n_starts", Int(alt.n_starts))
            setcolumn!(table, "phase_gibbs_J_per_mol_C", [Float64(get(alt.phase_gibbs, p, NaN)) for p in table["phase"]])
            setcolumn!(table, "start_labels", join(alt.start_labels, "; "))
            push!(tables, table)
        end
    end
    isempty(tables) && return DataTable(OrderedDict{String,Any}[]; columns=ALTERNATIVES_FIRST)
    table = concat(tables)
    return select(table, vcat(ALTERNATIVES_FIRST, [c for c in table.columns if !(c in ALTERNATIVES_FIRST)]))
end

"""
    save(sweep, directory) -> directory

Write summary.csv, species.csv, phases.csv, alternatives.csv and case.toml
(plus comparison.csv, if the case has comparison data) to `directory`, as the
Python app does. Not exported (`HTGEquilibrium.save`).
"""
function save(s::SweepResult, directory::AbstractString)
    mkpath(directory)
    to_csv(summary(s), joinpath(directory, "summary.csv"))
    to_csv(species_table(s), joinpath(directory, "species.csv"))
    to_csv(phase_summary(s), joinpath(directory, "phases.csv"))
    to_csv(alternatives_table(s), joinpath(directory, "alternatives.csv"))
    saved = replace(s.case)
    if py_truth(s.case.comparison) && isfile(String(s.case.comparison["file"]))
        cp(String(s.case.comparison["file"]), joinpath(directory, "comparison.csv"); force=true)
        comparison = OrderedDict{String,Any}(s.case.comparison)
        comparison["file"] = "comparison.csv"
        saved.comparison = comparison
    end
    save(saved, joinpath(directory, "case.toml"))
    return directory
end

# ---------------------------------------------------------------- chosen results
_point(T, P) = (py_float(T), py_float(P))

"""
    read_choices(folder) -> OrderedDict{(T_C, P_bar) => label}

The result each point uses, from choices.toml in `folder` (empty without the
file). The file holds only points that do not use result A:

    [[point]]
    T_C = 400
    P_bar = 240
    result = "B"
"""
function read_choices(folder)
    out = OrderedDict{Tuple{Float64,Float64},String}()
    path = joinpath(folder, CHOICES_FILE)
    isfile(path) || return out
    for entry in get(TOML.parsefile(path), "point", Any[])
        label = py_str(get(entry, "result", "A"))
        label == "A" || (out[_point(entry["T_C"], entry["P_bar"])] = label)
    end
    return out
end

"""
    write_choices(folder, choices) -> path or nothing

Write choices.toml (points that use A are left out; no file, and an old one
removed, when every point uses A).
"""
function write_choices(folder, choices)
    path = joinpath(folder, CHOICES_FILE)
    chosen = OrderedDict{Tuple{Float64,Float64},String}()
    for (point, label) in choices
        py_str(label) == "A" || (chosen[_point(point...)] = py_str(label))
    end
    if isempty(chosen)
        rm(path; force=true)
        return nothing
    end
    lines = ["# The result each point uses in plots and tables (points not listed use A).", ""]
    for (T, P) in sort!(collect(keys(chosen)))
        append!(lines, ["[[point]]", "T_C = $(py_repr(T))", "P_bar = $(py_repr(P))",
                        "result = \"$(chosen[(T, P)])\"", ""])
    end
    write(path, join(lines, "\n"))
    return path
end

_isclose_np(a, b) = abs(a - b) <= 1e-8 + 1e-5 * abs(b)   # numpy.isclose defaults
_tofloat(x) = x isa Real ? Float64(x) : x isa AbstractString ? something(py_tryfloat(x), NaN) : NaN

"""Rows of `table` at (T, P) (numpy isclose; P only when the table has P_bar)."""
function _at(table::DataTable, T, P)
    mask = [_isclose_np(_tofloat(x), T) for x in table["T_C"]]
    haskey(table, "P_bar") && (mask .&= [_isclose_np(_tofloat(x), P) for x in table["P_bar"]])
    return mask
end

_nansum(xs) = sum((Float64(x) for x in xs if !_isnull(x)); init=0.0)

function _reinfer!(t::DataTable, c)
    values, numeric = _infer_column(t.data[c])
    t.data[c] = values
    numeric ? push!(t.numeric, c) : delete!(t.numeric, c)
end

"""
    apply_choices(summary, species, phases, alternatives, choices) -> (summary, species, phases)

Tables with the chosen result at each point. `alternatives` is the
alternatives table (alternatives.csv); `choices` {(T, P) => label}. The summary
gains a "result" column when there are alternatives. At a point that uses
another result than A, its amounts replace the species rows, the phase rows
are recomputed from them, and the summary's Gibbs energy and gas carbon
follow; the solver columns stay those of A. Unknown labels are ignored.
"""
function apply_choices(summary::DataTable, species::DataTable, phases::DataTable, alternatives, choices)
    summary, species, phases = copy(summary), copy(species), copy(phases)
    (alternatives === nothing || isempty(alternatives)) && return summary, species, phases
    if !haskey(summary, "result")
        status = haskey(summary, "status") ? summary["status"] : fill("solved", nrows(summary))
        insert_column!(summary, min(3, length(summary.columns)), "result",
                       [isequal(s, "solved") ? "A" : "" for s in status])
    end
    species_columns = [c for c in alternatives.columns if !(c in ALTERNATIVE_COLUMNS) && haskey(species, c)]
    for (point, label) in choices
        T, P = _point(point...)
        picked = findall(_at(alternatives, T, P) .& [isequal(x, label) for x in alternatives["result"]])
        (isempty(picked) || label == "A") && continue
        mask = findall(_at(species, T, P))
        length(mask) == length(picked) || continue
        for c in species_columns
            for (k, i) in zip(mask, picked)
                species.data[c][k] = alternatives[c][i]
            end
            _reinfer!(species, c)
        end
        phase_of = alternatives["phase"]
        carbon = OrderedDict{Any,Float64}()   # groupby(phase).carbon_mol_per_kg_dry_feed.sum()
        gibbs = OrderedDict{Any,Float64}()    # groupby(phase).phase_gibbs_J_per_mol_C.first()
        for i in picked
            _isnull(phase_of[i]) && continue
            c = alternatives["carbon_mol_per_kg_dry_feed"][i]
            carbon[phase_of[i]] = get(carbon, phase_of[i], 0.0) + (_isnull(c) ? 0.0 : Float64(c))
            g = haskey(alternatives, "phase_gibbs_J_per_mol_C") ? alternatives["phase_gibbs_J_per_mol_C"][i] : NaN
            !haskey(gibbs, phase_of[i]) && !_isnull(g) && (gibbs[phase_of[i]] = Float64(g))
        end
        total_carbon = sum(values(carbon); init=0.0)
        for k in findall(_at(phases, T, P))
            phase = phases["phase"][k]
            part = [i for i in picked if isequal(phase_of[i], phase)]
            phases.data["mol_per_mol_C"][k] = _nansum(alternatives["mol_per_mol_C"][part])
            haskey(phases, "g_per_kg_dry_feed") &&
                (phases.data["g_per_kg_dry_feed"][k] = _nansum(alternatives["g_per_kg_dry_feed"][part]))
            phases.data["carbon_percent"][k] = total_carbon != 0 ? 100 * get(carbon, phase, 0.0) / total_carbon : 0.0
            haskey(phases, "gibbs_J_per_mol_C") && (phases.data["gibbs_J_per_mol_C"][k] = get(gibbs, phase, 0.0))
        end
        for c in ("mol_per_mol_C", "g_per_kg_dry_feed", "carbon_percent", "gibbs_J_per_mol_C")
            haskey(phases, c) && _reinfer!(phases, c)
        end
        for k in findall(_at(summary, T, P))
            summary.data["result"][k] = label
            summary.data["gibbs_J_per_mol_C"][k] = Float64(alternatives["gibbs_J_per_mol_C"][picked[1]])
            haskey(summary, "gas_carbon_percent") &&
                (summary.data["gas_carbon_percent"][k] =
                     total_carbon != 0 ? 100 * get(carbon, "gas", 0.0) / total_carbon : 0.0)
        end
        for c in ("result", "gibbs_J_per_mol_C", "gas_carbon_percent")
            haskey(summary, c) && _reinfer!(summary, c)
        end
    end
    return summary, species, phases
end

# ---------------------------------------------------------------- saved runs (Python: htg/gui/runs.py, data part)
"""
    RunTables

Everything the result and comparison views need from one run. `summary`,
`species` and `phases` use the chosen result at each point (`choices`, from
choices.toml; result A where none is chosen); `base` holds the tables as saved
(all A); `alternatives` is the alternatives table.
"""
struct RunTables
    label::String
    case::Union{Nothing,Case}
    summary::DataTable
    species::DataTable
    phases::DataTable
    comparison::Union{Nothing,DataTable}
    path::Union{Nothing,String}
    alternatives::DataTable
    choices::OrderedDict{Tuple{Float64,Float64},String}
    base::Tuple{DataTable,DataTable,DataTable}
end

"""
    RunTables(label, case, summary, species, phases; comparison, path, alternatives, choices)

Tables with `choices` applied (the given tables are kept as `base`; Python: `RunTables.build`).
"""
function RunTables(label, case, summary, species, phases; comparison=nothing, path=nothing, alternatives=nothing,
                   choices=nothing)
    alternatives = alternatives === nothing ? DataTable() : alternatives
    tables = RunTables(String(label), case, summary, species, phases, comparison,
                       path === nothing ? nothing : String(path), alternatives,
                       OrderedDict{Tuple{Float64,Float64},String}(), (summary, species, phases))
    return with_choices(tables, choices === nothing ? Dict() : choices)
end

"""These tables with another set of chosen results ({(T, P) => label})."""
function with_choices(t::RunTables, choices)
    known = OrderedDict{Tuple{Float64,Float64},String}(_point(point...) => String(label)
                                                       for (point, label) in choices if label != "A")
    summary, species, phases = apply_choices(t.base..., t.alternatives, known)
    return RunTables(t.label, t.case, summary, species, phases, t.comparison, t.path, t.alternatives, known, t.base)
end

"""(T_C, P_bar) of the points with more than one result, in table order."""
function flagged_points(t::RunTables)
    points = Tuple{Float64,Float64}[]
    isempty(t.alternatives) && return points
    for (T, P) in zip(t.alternatives["T_C"], t.alternatives["P_bar"])
        point = (Float64(T), Float64(P))
        point in points || push!(points, point)
    end
    return points
end

"""
    load_saved(folder, results_dir=nothing) -> RunTables

A saved run (a folder with summary.csv, species.csv, phases.csv and optionally
alternatives.csv, case.toml, comparison.csv, choices.toml), written by either version.
"""
function load_saved(folder, results_dir=nothing)
    folder = String(folder)
    file(name) = joinpath(folder, name)
    case = isfile(file("case.toml")) ? load_case(file("case.toml")) : nothing
    comparison = isfile(file("comparison.csv")) ? read_csv(file("comparison.csv")) : nothing
    label = results_dir === nothing ? basename(rstrip(folder, '/')) : relpath(folder, results_dir)
    summary, species, phases = (read_csv(file(name)) for name in ("summary.csv", "species.csv", "phases.csv"))
    if !isempty(phases) && !haskey(phases, "P_bar") && haskey(summary, "P_bar")
        # results saved before pressure grids existed (one pressure)
        insert_column!(phases, 1, "P_bar", summary["P_bar"][1])
    end
    alternatives = isfile(file("alternatives.csv")) ? read_csv(file("alternatives.csv")) : DataTable()
    return RunTables(label, case, summary, species, phases; comparison, path=folder, alternatives,
                     choices=read_choices(folder))
end
