# Results of equilibrium calculations (Python: htg/results.py, the data part).
#
# Amounts are in mol per mol of dry-feed carbon. The tables (species_table,
# phase_summary, summary, SweepResult, alternatives table, CSV files and
# choices.toml) are a later step; this file holds the data a solve fills.

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
