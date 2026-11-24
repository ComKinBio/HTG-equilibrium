# Constrained Gibbs-energy minimization at one temperature and pressure
# (Python: htg/equilibrium.py, and Model.solve of htg/model.py).
#
#     min_n  G(n) = sum over phases G_phase(n_phase)
#     s.t.   elemental balances, aqueous electroneutrality,
#            optional carbon-allocation constraints, n >= lower bounds
#
# Every (phase, species) pair is an independent variable in mol per mol of
# dry-feed carbon. A feasibility LP gives the starting point; the selected
# optimizer (default IPOPT) minimizes G. A converged result is a local minimum
# of the selected-phase problem; it is not certified to be the global one.
#
# Thread safety: `solve` builds everything it changes (phase objects, the
# GibbsProblem, the IPOPT problem) itself; the Model is only read. Several
# points can be solved at once on different threads.

"""The calculation failed; the message starts with the stage that failed (Python: SolveError(RuntimeError))."""
struct SolveError <: Exception
    msg::String
    cause::Any
end
SolveError(msg::AbstractString) = SolveError(String(msg), nothing)
Base.showerror(io::IO, e::SolveError) = print(io, e.msg)

"""
Python's class name of an error, for messages such as
"feed and initial composition: ValueError: ...". Julia types that stand for a
plain Python ValueError (CaseError, PhaseError, EquilibriumError, the data
layer's ArgumentErrors, ...) give "ValueError".
"""
py_exception_name(::EOSFailure) = "EOSFailure"
py_exception_name(::FeedError) = "FeedError"
py_exception_name(::LegacyConversionError) = "LegacyConversionError"
py_exception_name(::Union{ValueError,ArgumentError,FormulaError}) = "ValueError"
py_exception_name(::KeyError) = "KeyError"
py_exception_name(::RuntimeError) = "RuntimeError"
py_exception_name(::SolveError) = "SolveError"
py_exception_name(::DivideError) = "ZeroDivisionError"
py_exception_name(e) = string(nameof(typeof(e)))

_is_list(x) = x isa Union{AbstractVector,Tuple}

"""Elemental wt% of a dry feed given as relative molar amounts of registered species."""
function species_to_wt_percent(registry, species_moles, atomic_weights)
    entries = merge(all_molecular(registry), registry.solid, registry.ion)
    atoms = OrderedDict{String,Float64}()
    for (name, amount) in species_moles
        if !haskey(entries, name) || !isfinite(amount) || amount < 0
            throw(EquilibriumError("Invalid dry-feed species or amount: $(py_str(name))=$(py_str(amount))"))
        end
        for (element, count) in entries[name]["elements"]
            atoms[element] = get(atoms, element, 0.0) + amount * count
        end
    end
    masses = OrderedDict{String,Float64}(e => count * atomic_weights[e] for (e, count) in atoms)
    total = sum(values(masses); init=0.0)
    (total <= 0 || get(masses, "C", 0) <= 0) &&
        throw(EquilibriumError("Dry feed must contain positive mass and carbon"))
    return OrderedDict{String,Float64}(e => mass / total * 100 for (e, mass) in masses)
end

"""
    apply_phase_rules(case, registry) -> (phases, minimum, initial, liquid, gas_water, report)

Effective phases, bounds, start and initialization flags at the case temperature.

The aqueous cut-off is a modelling rule (not a computed phase boundary): at
T >= aqueous_cutoff_C, or with aqueous_enabled = false, the aqueous phase and
its bounds are removed. Feed atoms are kept and redistributed. The cut-off
takes priority over an aqueous carbon bound: the bound is dropped and
report["dropped_constraints"] says so (["aqueous carbon"]). `registry` needs
only `ion` (the ion table).
"""
function apply_phase_rules(case, registry)
    phases = OrderedDict{String,Vector{String}}(phase => collect(String, get(case.phases, phase, String[]))
                                                for phase in PHASES)
    disabled = case.aqueous_enabled === false
    cutoff = case.aqueous_cutoff_C
    if cutoff !== nothing
        cutoff = py_float(cutoff)
        (!isfinite(cutoff) || cutoff <= -273.15) &&
            throw(EquilibriumError("aqueous_cutoff_C must be finite and above absolute zero"))
        disabled = disabled || py_float(case.temperature_C) >= cutoff
    end
    liquid = case.initialization in ("auto", "liquid")
    gas_water = case.initialization == "vapor" || (case.initialization == "auto" && disabled)
    minimum = OrderedDict{Pair2,Float64}(case.minimum_amounts)
    removed = disabled ? copy(phases["aqueous"]) : String[]
    if disabled
        phases["aqueous"] = String[]
        filter!(kv -> kv[1][1] != "aqueous", minimum)
        liquid = false
    end
    _, dropped = effective_carbon_bounds(case.carbon_bounds, disabled)
    for (phase, names) in phases
        phase == "aqueous" && continue
        forbidden = [s for s in names if haskey(registry.ion, s)]
        isempty(forbidden) ||
            throw(EquilibriumError("Ions are only allowed in aqueous: $phase: $(py_repr(forbidden))"))
    end
    report = OrderedDict{String,Any}("aqueous_disabled" => disabled, "aqueous_cutoff_C" => cutoff,
                                     "removed_aqueous_species" => removed, "initial_composition_reset" => false,
                                     "dropped_constraints" => dropped)
    initial = case.initial_amounts
    if initial !== nothing
        any(v -> !isfinite(v) || v < 0, values(initial)) && throw(EquilibriumError("Invalid initial amount."))
        selected = Set(Pair2[(p, s) for (p, names) in phases for s in names])
        lost = [pair for (pair, v) in initial if !(pair in selected) && v > 0]
        if !isempty(lost)
            initial = nothing
            report["initial_composition_reset"] = true
            report["initial_reset_reason"] = "Populated species removed: $(py_repr(lost))"
        else
            initial = OrderedDict{Pair2,Float64}(pair => v for (pair, v) in initial if pair in selected)
        end
    end
    return phases, minimum, initial, liquid, gas_water, report
end

"""Name of the calculation stage, reported when something fails."""
mutable struct Stage
    name::String
end
Stage() = Stage("case definition")

"""
    Setup

Everything needed to evaluate G and run the optimizer for one case (Python:
the `Setup` dataclass; same fields). Positions (`free`, `indices`) are 1-based;
`phase_G(amounts)` and `total_G(amounts)` evaluate the phase models.
"""
struct Setup
    case::Case
    pairs::Vector{Pair2}
    phase_species::OrderedDict{String,Vector{String}}
    indices::OrderedDict{String,Vector{Int}}
    phase_models::OrderedDict{String,Any}
    phase_G::Function
    total_G::Function
    initial_amounts::Vector{Float64}
    initial_G::Float64
    x_start::Vector{Float64}
    free::Vector{Int}
    bounds::Vector{Tuple{Float64,Nothing}}
    A_opt::Matrix{Float64}
    b_opt::Vector{Float64}
    A_full::Matrix{Float64}
    b_full::Vector{Float64}
    A::Matrix{Float64}
    b::Vector{Float64}
    charge::Vector{Float64}
    gas_carbon::Vector{Float64}
    feed_carbon::Float64
    ineq_rows::Matrix{Float64}      # rows . x >= targets on the free amounts (min and -max bounds)
    ineq_targets::Vector{Float64}
    carbon_checks::Vector{Tuple{String,String,Vector{Float64},Float64}}  # (phase, kind, full row, fraction)
    energy_scale::Float64
    elements::Vector{String}
    entries::Vector{Any}
    C_mol_per_g::Float64
    T::Float64
    atomic_weights::Dict{String,Float64}
    phase_rules::OrderedDict{String,Any}
end

function expand(s::Setup, x)
    amounts = zeros(length(s.pairs))
    amounts[s.free] .= x
    return amounts
end

"""
    setup(model, case, stage=Stage(); build_phases=true) -> Setup

Constraints, feasible start and phase models for `case` (stops before
optimizing). With build_phases=false it stops after the feasibility LP (no
equation of state) and returns `nothing`. The electrolyte phase uses the
case's water model (`water_model_of(case)`).
"""
function setup(model::Model, case::Case, stage::Stage=Stage(); build_phases=true)
    validate(case)
    if case.points !== nothing || _is_list(case.temperature_C) || _is_list(case.pressure_bar)
        throw(EquilibriumError("A single temperature and pressure is needed; use Model.sweep for a grid, " *
                               "or case.at(T, P)."))
    end
    registry = model.registry
    stage.name = "feed and initial composition"
    AW = Dict{String,Float64}(ATOMIC_WEIGHTS)
    for (e, w) in case.atomic_weights
        AW[String(e)] = Float64(w)
    end
    phase_species, minimum_amounts, initial_input, liquid_init, gas_water_init, allocation =
        apply_phase_rules(case, registry)
    TEMPERATURE_C = py_float(case.temperature_C)
    PRESSURE_BAR = py_float(case.pressure_bar)
    feed = resolved(case.feed, molar_mass_function(registry_elements(registry), AW))
    FEED = dry_element_wt_percent(feed, moles -> species_to_wt_percent(registry, moles, AW))
    WATER = total_water_g_per_g_dry(feed)
    SALTS = feed.salts
    BOUNDS, _ = effective_carbon_bounds(case.carbon_bounds, allocation["aqueous_disabled"])
    WATER_FLOOR = py_float(case.aqueous_water_floor)
    AQUEOUS_MODEL = case.aqueous_model
    if AQUEOUS_MODEL == "molecular" && any(s -> haskey(registry.ion, s), phase_species["aqueous"])
        throw(EquilibriumError("Ions require the electrolyte aqueous model"))
    end
    INITIAL_GAS_WATER_FRACTION = py_float(case.initial_gas_water_fraction)
    (!isfinite(INITIAL_GAS_WATER_FRACTION) || INITIAL_GAS_WATER_FRACTION < 0) &&
        throw(EquilibriumError("Invalid initial gas-water fraction"))
    (!isfinite(TEMPERATURE_C) || TEMPERATURE_C <= -273.15) && throw(EquilibriumError("Invalid temperature."))
    (!isfinite(PRESSURE_BAR) || PRESSURE_BAR <= 0) && throw(EquilibriumError("Pressure must be finite and positive."))
    (!isfinite(WATER) || WATER < 0) && throw(EquilibriumError("Added water must be finite and nonnegative."))
    any(v -> !isfinite(v) || v < 0, values(FEED)) && throw(EquilibriumError("Invalid feed elemental composition."))
    get(FEED, "C", 0.0) <= 0 &&
        throw(EquilibriumError("This carbon-normalized calculation requires positive dry-feed carbon."))
    sum(values(FEED); init=0.0) > 100.0 + 1e-08 && throw(EquilibriumError("Feed elemental percentages exceed 100%."))
    if !isempty(phase_species["aqueous"]) && AQUEOUS_MODEL == "electrolyte" &&
       (!isfinite(WATER_FLOOR) || WATER_FLOOR <= 1e-09)
        throw(EquilibriumError("aqueous_water_floor must exceed the aqueous model's 1e-9 mol solvent cutoff."))
    end
    for (phase, names) in phase_species
        length(names) != length(unique(names)) && throw(EquilibriumError("Duplicate species in $phase."))
    end

    pairs = Pair2[(phase, s) for (phase, names) in phase_species for s in names]
    isempty(pairs) && throw(EquilibriumError("Select at least one species."))
    entries = Any[entry(registry, phase, s) for (phase, s) in pairs]
    indices = OrderedDict{String,Vector{Int}}(phase => [j for (j, (p, _)) in enumerate(pairs) if p == phase]
                                              for phase in PHASES)
    element_set = Set{String}(keys(FEED))
    push!(element_set, "H", "O")
    for salt in values(SALTS), e in keys(salt.elements)
        push!(element_set, String(e))
    end
    for ent in entries, e in keys(ent["elements"])
        push!(element_set, String(e))
    end
    elements = sort!(collect(element_set))
    missing_aw = [e for e in elements if !haskey(AW, e)]
    isempty(missing_aw) || throw(EquilibriumError("Missing atomic weights: " * join(missing_aw, ", ")))
    A = Float64[Float64(get(ent["elements"], e, 0.0)) for e in elements, ent in entries]
    C_mol_per_g = FEED["C"] / 100.0 / AW["C"]
    b = Float64[get(FEED, e, 0.0) / 100.0 / AW[e] / C_mol_per_g for e in elements]
    water_MW = 2.0 * AW["H"] + AW["O"]
    water_mol_per_C = WATER / water_MW / C_mol_per_g
    iH, iO = findfirst(==("H"), elements), findfirst(==("O"), elements)
    b[iH] += 2.0 * water_mol_per_C
    b[iO] += water_mol_per_C
    for (salt_name, salt) in SALTS
        amount = py_float(salt.mol_per_g_dry_feed)
        (!isfinite(amount) || amount < 0) && throw(EquilibriumError("Invalid amount for $salt_name."))
        for (element, count) in salt.elements
            (!isfinite(count) || count < 0) && throw(EquilibriumError("Invalid elemental count for $salt_name."))
            b[findfirst(==(String(element)), elements)] += amount * count / C_mol_per_g
        end
    end
    iC = findfirst(==("C"), elements)
    carbon = A[iC, :]
    feed_carbon = b[iC]
    gas_carbon = carbon .* Float64[phase == "gas" for (phase, _) in pairs]
    charge = Float64[phase == "aqueous" && haskey(registry.ion, s) ? Float64(registry.ion[s]["charge"]) : 0.0
                     for (phase, s) in pairs]
    missing_ = [e for (k, e) in enumerate(elements) if b[k] > 0 && !any(>(0), A[k, :])]
    isempty(missing_) || throw(EquilibriumError("No allowed product carrier for feed elements: " *
                                                join(missing_, ", ") *
                                                ". Add neutral gas species or solids that contain them."))
    # Carbon bounds: exact -> equality rows; min / max -> inequality rows.
    exact_rows, bound_rows = carbon_bound_rows(BOUNDS, pairs, carbon, case.solid_carbon_target)
    exact = OrderedDict{String,Float64}(phase => fraction for (phase, _, fraction) in exact_rows)
    if haskey(exact, "gas") && haskey(exact, "solid") && exact["gas"] + exact["solid"] > 1.0 + 1e-10
        throw(EquilibriumError("Gas-carbon plus solid-carbon fractions exceed one."))
    end
    if sum(values(exact); init=0.0) + sum((f for (_, kind, _, f) in bound_rows if kind == "min"); init=0.0) >
       1.0 + 1e-10
        throw(EquilibriumError("The carbon shares add up to more than 100 %."))
    end
    rows = Vector{Float64}[A[k, :] for k in axes(A, 1)]
    push!(rows, charge)
    targets = vcat(b, 0.0)
    carbon_constraints = Tuple{Vector{Float64},Float64}[]
    carbon_checks = Tuple{String,String,Vector{Float64},Float64}[]
    for (phase, row, fraction) in exact_rows
        push!(rows, row)
        push!(targets, fraction * feed_carbon)
        push!(carbon_constraints, (row, fraction))
        push!(carbon_checks, (phase, "exact", row, fraction))
    end
    inequalities = Tuple{String,String,Vector{Float64},Float64}[]
    for (phase, kind, row, fraction) in bound_rows
        push!(carbon_checks, (phase, kind, row, fraction))
        if kind == "max" && fraction == 0.0
            push!(carbon_constraints, (row, 0.0))  # no carbon there: those species are fixed at zero
            continue
        end
        push!(inequalities, (phase, kind, row, fraction))
    end
    A_full = permutedims(reduce(hcat, rows))
    b_full = targets
    lower = zeros(length(pairs))
    for (pair, value) in minimum_amounts
        j = findfirst(==(pair), pairs)
        j === nothing && throw(EquilibriumError("Lower-bound species not selected: $(py_repr(pair))"))
        (!isfinite(value) || value < 0) && throw(EquilibriumError("Invalid lower bound for $(py_repr(pair))."))
        lower[j] = value
    end
    if !isempty(phase_species["aqueous"])
        "water" in phase_species["aqueous"] || throw(EquilibriumError("The aqueous phase requires solvent water."))
        w = findfirst(==(("aqueous", "water")), pairs)
        lower[w] = max(lower[w], WATER_FLOOR)
    end
    fixed_zero = falses(length(pairs))
    # Species containing an element absent from the feed must be zero.
    for k in axes(A, 1)
        b[k] == 0.0 && (fixed_zero .|= A[k, :] .> 0)
    end
    for (row, fraction) in carbon_constraints
        if fraction == 0.0
            fixed_zero .|= row .> 0
        elseif fraction == 1.0
            fixed_zero .|= (carbon .> 0) .& (row .== 0)
        end
    end
    any(lower[fixed_zero] .> 0) && throw(EquilibriumError("A minimum amount conflicts with a carbon constraint."))
    free = findall(.!fixed_zero)
    isempty(free) && throw(EquilibriumError("No free species remain."))
    expand_(x) = (amounts = zeros(length(pairs)); amounts[free] .= x; amounts)

    reduced = A_full[:, free]
    F = qr(permutedims(reduced), ColumnNorm())   # LAPACK geqp3, as scipy.linalg.qr(pivoting=True)
    diagonal = abs.(diag(F.R))
    rank = count(diagonal .> 1e-10 * max(1.0, isempty(diagonal) ? 0.0 : maximum(diagonal)))
    selected_rows = sort(F.p[1:rank])
    A_opt = reduced[selected_rows, :]
    b_opt = b_full[selected_rows]
    bounds = Tuple{Float64,Nothing}[(lower[j], nothing) for j in free]
    start_objective = zeros(length(free))
    if gas_water_init
        for (k, j) in enumerate(free)
            pairs[j] == ("gas", "water") && (start_objective[k] = -1.0)
        end
    end
    start_upper = fill(Inf, length(free))
    if liquid_init && ("aqueous", "water") in pairs
        for (k, j) in enumerate(free)
            if pairs[j] == ("aqueous", "water")
                start_objective[k] = -1.0
            elseif pairs[j] == ("gas", "water")
                start_upper[k] = max(lower[j], INITIAL_GAS_WATER_FRACTION * water_mol_per_C)
            end
        end
    end
    # min and max bounds are inequalities in the LP, the optimizer and the final check,
    # written as rows . x >= targets (a max bound as -row . x >= -max).
    signs = [kind == "min" ? 1.0 : -1.0 for (_, kind, _, _) in inequalities]
    ineq_full = zeros(length(inequalities), length(pairs))
    for (k, (sign, (_, _, row, _))) in enumerate(zip(signs, inequalities))
        ineq_full[k, :] .= sign .* row
    end
    ineq_rows = ineq_full[:, free]
    ineq_targets = Float64[sign * fraction * feed_carbon for (sign, (_, _, _, fraction)) in zip(signs, inequalities)]
    has_ineq = !isempty(inequalities)
    success, x_lp, message = linprog(start_objective; A_eq=A_opt, b_eq=b_opt, lower=lower[free], upper=start_upper,
                                     A_ub=has_ineq ? -ineq_rows : nothing, b_ub=has_ineq ? -ineq_targets : nothing)
    success || throw(EquilibriumError("Selected feed, species and constraints are infeasible: " * message))
    x_start = x_lp
    initial_amounts = expand_(x_lp)
    if initial_input !== nothing
        initial_amounts = Float64[get(initial_input, pair, 0.0) for pair in pairs]
        isempty(setdiff(keys(initial_input), pairs)) ||
            throw(EquilibriumError("Initial composition contains an unselected phase/species"))
        if any(!isfinite, initial_amounts) || any(initial_amounts .< lower) || any(initial_amounts[fixed_zero] .!= 0)
            throw(EquilibriumError("Initial composition violates species bounds"))
        end
        x_start = initial_amounts[free]
    end
    if maximum(abs.(A_full * initial_amounts .- b_full)) > 1e-08
        throw(EquilibriumError("The selected species cannot satisfy all original balances and carbon constraints."))
    end
    for (k, (phase, kind, _, _)) in enumerate(inequalities)
        if ineq_full[k, :]' * initial_amounts < ineq_targets[k] - 1e-8
            label = lowercase(PHASE_LABELS[phase])
            throw(EquilibriumError("Initial composition violates $(kind)imum $label-carbon fraction"))
        end
    end

    build_phases || return nothing
    stage.name = "thermochemistry / aqueous validity"
    T = TEMPERATURE_C + 273.15
    P = PRESSURE_BAR * 100000.0
    electrolyte = !isempty(phase_species["aqueous"]) && AQUEOUS_MODEL == "electrolyte"
    to_resolve = String[s for (phase, names) in phase_species if phase != "solid" for s in names
                        if !haskey(registry.ion, s)]
    electrolyte && append!(to_resolve, registry.aqueous_reference_neutrals)
    params = resolve_parameters(registry, model.databases, T, unique(to_resolve); atomic_weights=AW)
    built = OrderedDict{String,Any}()
    if electrolyte
        built["aqueous"] = ElectrolytePhase(water_model_of(case), model, phase_species["aqueous"], T, P, params, AW)
    end
    for (phase, names) in phase_species
        (isempty(names) || haskey(built, phase)) && continue
        if phase == "solid"
            built[phase] = SolidPhase(names, solid_g0(model, names, T, P))
        else
            root = phase == "gas" ? "vapour" : "liquid"
            built[phase] = FluidPhase(phase, names, nasa_g0(model, names, T),
                                      build_pcsaft_model(names, params, registry.binary_interactions), T, P, root)
        end
    end
    phase_models = OrderedDict{String,Any}(phase => built[phase] for phase in PHASES if haskey(built, phase))

    function phase_G(amounts)
        terms = OrderedDict{String,Float64}()
        for phase in PHASES
            n = amounts[indices[phase]]
            terms[phase] = isempty(phase_species[phase]) || sum(n; init=0.0) <= 0.0 ? 0.0 :
                           Float64(gibbs(phase_models[phase], n))
        end
        all(isfinite, values(terms)) || throw(EquilibriumError("Nonfinite phase Gibbs energy."))
        return terms
    end
    function total_G(amounts)
        minimum(amounts) < -1e-10 && throw(EquilibriumError("Negative species amount."))
        return sum(values(phase_G(max.(amounts, 0.0))))
    end

    stage.name = "initial Gibbs evaluation"
    initial_G = total_G(initial_amounts)
    energy_scale = case.energy_scale !== nothing ? py_float(case.energy_scale) : max(10000.0, abs(initial_G))
    (!isfinite(energy_scale) || energy_scale <= 0) && throw(EquilibriumError("energy_scale must be finite and positive"))
    return Setup(case, pairs, phase_species, indices, phase_models, phase_G, total_G, initial_amounts, initial_G,
                 x_start, free, bounds, A_opt, b_opt, A_full, b_full, A, b, charge, gas_carbon, feed_carbon,
                 ineq_rows, ineq_targets, carbon_checks, energy_scale, elements, entries, C_mol_per_g, T, AW,
                 allocation)
end

"""
    carbon_bound_rows(bounds, pairs, carbon, solid_carbon_target) -> (exact_rows, bound_rows)

Linear rows on phase carbon for the carbon bounds.

`bounds` is {phase => {"exact" => f} or {"min" => f, "max" => f}} (fractions
of feed carbon); `pairs` [(phase, species)]; `carbon` the carbon atoms of each
pair. Returns [(phase, row, fraction)] for equalities (row . n = fraction *
feed carbon) and [(phase, kind, row, fraction)] for "min" (row . n >= ...)
and "max" (row . n <= ...), in the order of `bounds`. "solid" counts solid
carbon (char) only, or all solids (solid_carbon_target).
"""
function carbon_bound_rows(bounds, pairs, carbon, solid_carbon_target)
    carbon = collect(Float64, carbon)
    exact_rows = Tuple{String,Vector{Float64},Float64}[]
    bound_rows = Tuple{String,String,Vector{Float64},Float64}[]
    for (phase, bound) in bounds
        if phase == "solid" && solid_carbon_target == "char"
            j = findfirst(==(("solid", "carbon")), pairs)
            j === nothing && throw(EquilibriumError("Select solid 'carbon' to constrain char carbon."))
            row = zeros(length(pairs))
            row[j] = carbon[j]
        else
            row = carbon .* Float64[p == phase for (p, _) in pairs]
        end
        for kind in BOUND_KINDS
            haskey(bound, kind) || continue
            fraction = py_float(bound[kind])
            (!isfinite(fraction) || !(0 <= fraction <= 1)) &&
                throw(EquilibriumError("$(get(PHASE_LABELS, phase, phase)) carbon $kind must be in [0, 1]."))
            if kind == "exact"
                push!(exact_rows, (String(phase), row, fraction))
            else
                push!(bound_rows, (String(phase), kind, row, fraction))
            end
        end
    end
    return exact_rows, bound_rows
end

"""Largest violation of the carbon bounds, in fractions of the feed carbon."""
function carbon_bound_violation(checks, amounts, feed_carbon)
    worst = 0.0
    for (_, kind, row, fraction) in checks
        share = (row' * amounts) / feed_carbon
        excess = kind == "min" ? fraction - share : kind == "max" ? share - fraction : abs(share - fraction)
        worst = max(worst, excess)
    end
    return worst
end

function _optimize(model::Model, s::Setup)
    case, registry = s.case, model.registry
    # Python's slsqp-fd (finite-difference SLSQP) and trust-constr are not ported.
    case.solver in UNAVAILABLE_SOLVERS && throw(EquilibriumError(unavailable_solver_message(case.solver)))
    pairs = s.pairs
    eos_errors = String[]
    has_ineq = size(s.ineq_rows, 1) > 0
    problem = GibbsProblem(
        phases=s.phase_models, pairs=pairs, free=s.free, phase_indices=s.indices, ions=Set(keys(registry.ion)),
        T=s.T, A_eq=s.A_opt, b_eq=s.b_opt, ineq_rows=has_ineq ? s.ineq_rows : nothing,
        ineq_targets=has_ineq ? s.ineq_targets : nothing, lower=[low for (low, _) in s.bounds],
        x_start=s.x_start, G0=s.initial_G, scale=s.energy_scale, model_total_G=s.total_G)
    options = OrderedDict{String,Any}(
        "method" => case.solver, "max_iterations" => py_int(case.max_iterations), "starts" => py_int(case.starts),
        "seed" => py_int(case.seed), "distinct_composition_tol" => py_float(case.distinct_composition_tol),
        "distinct_gibbs_tol" => py_float(case.distinct_gibbs_tol))
    for (k, v) in case.solver_options
        options[String(k)] = v
    end
    case.stability_test === nothing || (options["stability_test"] = py_truth(case.stability_test))
    optimized, diagnostics = optimize(problem, options)
    groups = haskey(diagnostics, "groups") ? pop!(diagnostics, "groups") : OrderedDict{String,Any}[]
    diagnostics["results"] = [OrderedDict{String,Any}(key => g[key] for key in ("label", "G", "n_starts", "start_labels"))
                              for g in groups]
    append!(eos_errors, problem.eos_errors[1:min(end, max(0, 10 - length(eos_errors)))])
    amounts = expand(s, optimized.x)
    residual = maximum(abs.(s.A_full * amounts .- s.b_full))
    violation = carbon_bound_violation(s.carbon_checks, amounts, s.feed_carbon)
    for check in s.carbon_checks
        phase, kind = check[1], check[2]
        kind == "exact" && continue  # equalities are checked with the balances
        excess = carbon_bound_violation([check], amounts, s.feed_carbon) * s.feed_carbon
        if excess > 1e-7
            label = lowercase(PHASE_LABELS[phase])
            throw(RuntimeError("$(kind == "min" ? "Minimum" : "Maximum") $label-carbon allocation " *
                               "violated by $(@sprintf("%.3e", excess)) mol C/mol dry-feed C"))
        end
    end
    if !optimized.success || any(!isfinite, amounts) || minimum(amounts) < -1e-09 || residual > 1e-07
        throw(RuntimeError("$(optimized.message); constraint residual=$(@sprintf("%.3e", residual)); " *
                           "recorded EOS errors=$(py_repr(eos_errors))"))
    end
    amounts = max.(amounts, 0.0)
    G_terms = s.phase_G(amounts)
    final_G = sum(values(G_terms))
    final_G > s.initial_G + 0.0001 && throw(RuntimeError("Final Gibbs energy exceeds the feasible start."))
    alternatives = Alternative[]
    for group in groups
        if group["label"] == "A"
            alt_amounts, terms, gibbs_ = amounts, G_terms, final_G
        else
            alt_amounts = max.(expand(s, group["x"]), 0.0)
            terms = try
                s.phase_G(alt_amounts)
            catch exc
                exc isa Union{ValueError,ArgumentError} || rethrow()
                continue
            end
            gibbs_ = sum(values(terms))
        end
        push!(alternatives, Alternative(group["label"], gibbs_, OrderedDict(zip(pairs, alt_amounts)),
                                        group["n_starts"], group["start_labels"]; phase_gibbs=terms))
    end
    species_data = OrderedDict{Pair2,OrderedDict{String,Any}}()
    for (pair, ent) in zip(pairs, s.entries)
        species_data[pair] = OrderedDict{String,Any}(
            "formula" => get(ent, "formula", pair[2]), "elements" => OrderedDict{String,Any}(ent["elements"]),
            "charge" => haskey(registry.ion, pair[2]) ? get(ent, "charge", 0) : 0)
    end
    return Result(
        temperature_C=py_float(case.temperature_C), pressure_bar=py_float(case.pressure_bar),
        amounts=OrderedDict{Pair2,Float64}(zip(pairs, amounts)),
        initial_amounts=OrderedDict{Pair2,Float64}(zip(pairs, s.initial_amounts)),
        gibbs_J_per_mol_C=final_G, phase_gibbs_J_per_mol_C=G_terms,
        balance_residual=residual, charge_residual=s.charge' * amounts,
        iterations=optimized.nit, eos_trial_errors=eos_errors,
        feed_C_mol_per_g_dry=s.C_mol_per_g, elements=s.elements, feed_atoms=copy(s.b),
        product_atoms=s.A * amounts, gas_carbon_fraction=(s.gas_carbon' * amounts) / s.feed_carbon,
        carbon_bound_violation=violation, phase_rules=s.phase_rules,
        optimizer=diagnostics, atomic_weights=s.atomic_weights, species_data=species_data,
        alternatives=alternatives, water_model=String(case.water_model))
end

"""Solve a single-point `case` (Python: htg.equilibrium.solve); failures are SolveErrors."""
function _solve_point(model::Model, case::Case)
    stage = Stage()
    try
        s = setup(model, case, stage)
        stage.name = "optimization / final validation"
        return _optimize(model, s)
    catch exc
        (exc isa SolveError || exc isa InterruptException) && rethrow()
        throw(SolveError("$(stage.name): $(py_exception_name(exc)): $(_exception_text(exc))", exc))
    end
end

"""
    solve(model, case, temperature_C=nothing, pressure_bar=nothing) -> Result

Equilibrium for `case` at one point (Python: `Model.solve`). Give
temperature_C (and pressure_bar, if the case has several pressures) to pick a
point of a sweep case; a single-point case needs neither. A failed calculation
throws a `SolveError` whose message starts with the stage that failed.

Thread safe: several points (or cases) may be solved at once on different
threads with one shared `Model`.
"""
function solve(model::Model, case::Case, temperature_C=nothing, pressure_bar=nothing)
    if temperature_C !== nothing || pressure_bar !== nothing
        if temperature_C === nothing
            temps = temperatures(case)
            length(temps) != 1 &&
                throw(EquilibriumError("This case has several temperatures; give temperature_C too."))
            temperature_C = temps[1]
        end
        case = at(case, temperature_C, pressure_bar)
    elseif case.points !== nothing || _is_list(case.temperature_C) || _is_list(case.pressure_bar)
        points = grid(case)
        length(points) != 1 &&
            throw(EquilibriumError("This case has several points; give temperature_C (and pressure_bar), " *
                                   "or use Model.sweep."))
        case = at(case, points[1]...)
    end
    return _solve_point(model, case)
end
