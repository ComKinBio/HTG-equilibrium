# Shared setup of the solver tests (test_optimizer.jl, test_equilibrium.jl):
# tests/reference/data/solves.json and the port of tests/reference/convert.py
# (reference entries in the old settings format -> Case). Uses Main.DataTestSetup
# (the frozen reference models). Included once into Main as `SolvesTestSetup`.
module SolvesTestSetup

using HTGEquilibrium
using HTGEquilibrium: OrderedDict
using Main.DataTestSetup: MODELS, reference

export SOLVES, SOLVED, FAILED, REFERENCE, LEGACY_PHASE, case_from_legacy, legacy_pairs, legacy_terms,
       FULL_TESTS

"""HTG_FULL_TESTS=1 runs every reference solve (and the slow checks); default: a representative subset."""
const FULL_TESTS = get(ENV, "HTG_FULL_TESTS", "") in ("1", "true", "yes")

const SOLVES = reference("solves")["cases"]
const SOLVED = [c for c in SOLVES if haskey(c, "result")]
const FAILED = [c for c in SOLVES if haskey(c, "error")]
const REFERENCE = Dict(c["name"] => c for c in SOLVES)

const PHASE = Dict("gas" => "gas", "aqueous" => "aqueous", "organic liquid" => "organic_liquid", "solid" => "solid")
const LEGACY_PHASE = Dict(v => k for (k, v) in PHASE)

"""[[legacy phase, species, value], ...] -> {(phase, species) => value} (convert.py `pairs`)."""
legacy_pairs(triples) = OrderedDict{Tuple{String,String},Float64}(
    (PHASE[p], String(s)) => Float64(v) for (p, s, v) in something(triples, Any[]))

"""Phase-Gibbs terms with the reference's (legacy) phase names."""
legacy_terms(terms) = Dict(LEGACY_PHASE[p] => v for (p, v) in terms)

_get(d, k, default) = (v = get(d, k, default); v === nothing ? default : v)

"""convert.py `case_from_legacy`: the reference entry as a Case (solver "slsqp-fd" when the
entry has no optimizer settings, as in Python)."""
function case_from_legacy(entry)
    s = entry["settings"]
    phases = OrderedDict{String,Any}(PHASE[p] => collect(String, names) for (p, names) in s["phase_species"])
    for (phase, names) in _get(s, "ion_product_species", Dict())
        phases[PHASE[phase]] = unique(vcat(get(phases, PHASE[phase], String[]), collect(String, names)))
    end
    flags = (Bool(_get(s, "liquid_initialization", false)), Bool(_get(s, "gas_water_initialization", false)))
    initialization = Dict((true, false) => "liquid", (false, true) => "vapor", (false, false) => "lp")[flags]
    optimizer = get(s, "optimizer", nothing)
    species_moles = get(s, "feed_species_moles", nothing)
    feed = Feed(elements_wt_percent=OrderedDict{String,Any}(_get(s, "feed_element_wt_percent", Dict())),
                species_moles=species_moles === nothing ? nothing : OrderedDict{String,Any}(species_moles),
                water_g_per_g_dry=Float64(_get(s, "water_g_per_g_dry_feed", 0.0)),
                salts=OrderedDict{String,Any}(k => OrderedDict{String,Any}(v) for (k, v) in _get(s, "added_salts", Dict())))
    initial = get(s, "initial_phase_amounts", nothing)
    return Case(
        name=entry["name"], temperature_C=Float64(s["temperature_C"]), pressure_bar=Float64(s["pressure_bar"]),
        feed=feed, phases=phases, atomic_weights=OrderedDict{String,Any}(entry["atomic_weights"]),
        aqueous_model=_get(s, "aqueous_model", "electrolyte"),
        aqueous_enabled=_get(s, "aqueous_enabled", true), aqueous_cutoff_C=get(s, "aqueous_cutoff_C", nothing),
        aqueous_water_floor=Float64(_get(s, "aqueous_water_floor", 1e-8)),
        gas_carbon_fraction=get(s, "gas_carbon_fraction", nothing),
        gas_carbon_min_fraction=get(s, "gas_carbon_min_fraction", nothing),
        solid_carbon_fraction=get(s, "solid_carbon_fraction", nothing),
        solid_carbon_target=Dict("CHAR" => "char", "all solids" => "all_solids")[_get(s, "solid_carbon_target", "CHAR")],
        phase_carbon_fractions=OrderedDict{String,Any}(PHASE[k] => v for (k, v) in _get(s, "phase_carbon_fractions", Dict())),
        minimum_amounts=legacy_pairs(get(s, "minimum_amounts", nothing)),
        initial_amounts=initial === nothing || isempty(initial) ? nothing : legacy_pairs(initial),
        initialization=initialization,
        initial_gas_water_fraction=Float64(_get(s, "initial_gas_water_fraction", 1e-6)),
        solver=optimizer === nothing ? "slsqp-fd" : optimizer["method"],
        multistart=optimizer === nothing ? 0 : Int(_get(optimizer, "multistart", 0)),
        max_iterations=Int(_get(s, "max_iterations", 5000)),
        fd_ftol=Float64(_get(s, "optimizer_ftol", 1e-9)), fd_eps=Float64(_get(s, "optimizer_eps", 1e-7)),
        energy_scale=get(s, "optimizer_energy_scale", nothing))
end

end # module
