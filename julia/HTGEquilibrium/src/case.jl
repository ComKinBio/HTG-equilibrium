# Case definition: everything that specifies one equilibrium calculation or
# sweep (Python: htg/case.py).
#
#     case = Case(temperature_C=450, pressure_bar=240,
#                 feed=Feed(elements_wt_percent=OrderedDict("C" => 30.4, ...), water_g_per_g_dry=4.0),
#                 phases=OrderedDict("gas" => [...], "aqueous" => [...], "solid" => ["carbon"]))
#     case = load_case("cases/manure_gas_min.toml")
#     grid(case), chains(case), temperatures(case), at(case, 300.0), replace(case; starts=3)
#
# Units: temperatures in °C, pressures in bar (absolute). Species amounts in
# constraints and results are mol per mol of dry-feed carbon.
#
# As in Python, the fields hold what the file gave (Int or Float64 numbers,
# OrderedDicts in file order) until `_post_init!` normalises the ones the
# calculation reads. Errors that are ValueErrors in Python are `CaseError`s
# (or `FeedError`s, `LegacyConversionError`s) with the same messages.

"""A case that cannot be used; the message says why (Python: ValueError)."""
struct CaseError <: ValueError
    msg::String
end

"""An old case file cannot be written in format 2; the message says why."""
struct LegacyConversionError <: ValueError
    msg::String
end

const SOLVERS = ("ipopt", "slsqp", "trust-constr", "slsqp-fd")
const INITIALIZATIONS = ("auto", "liquid", "vapor", "lp")

const _Dict = OrderedDict{String,Any}

# ---------------------------------------------------------------- Salt and Feed
"""An added salt: supplies atoms only; its identity does not constrain the products."""
mutable struct Salt
    mol_per_g_dry_feed::Any
    elements::Any
end
Salt(; mol_per_g_dry_feed, elements) = Salt(mol_per_g_dry_feed, elements)
Base.:(==)(a::Salt, b::Salt) = a.mol_per_g_dry_feed == b.mol_per_g_dry_feed && a.elements == b.elements

# Python's TypeErrors for Cls(**kwargs) with unknown or missing keys (Julia: ArgumentError, same text).
_unexpected(cls, key) = ArgumentError("$cls.__init__() got an unexpected keyword argument $(py_repr(string(key)))")
function _missing_arguments(cls, names)
    quoted = py_repr.(names)
    list = length(quoted) == 1 ? quoted[1] : length(quoted) == 2 ? join(quoted, " and ") :
           join(quoted[1:end-1], ", ") * ", and " * quoted[end]
    plural = length(names) == 1 ? "argument" : "arguments"
    return ArgumentError("$cls.__init__() missing $(length(names)) required positional $plural: $list")
end

_salt(s::Salt) = s
function _salt(d::AbstractDict)
    for k in keys(d)
        string(k) in ("mol_per_g_dry_feed", "elements") || throw(_unexpected("Salt", k))
    end
    missing_ = [k for k in ("mol_per_g_dry_feed", "elements") if !haskey(d, k)]
    isempty(missing_) || throw(_missing_arguments("Salt", missing_))
    return Salt(d["mol_per_g_dry_feed"], d["elements"])
end
_salts(salts) = OrderedDict{String,Salt}(string(name) => _salt(s) for (name, s) in salts)

"""
    Feed(; elements_wt_percent, basis, ash_wt_percent, unmodelled_wt_percent,
         species_moles, moisture_wt_percent, water_g_per_g_dry, salts, species_amounts)

Dry feed, water and salts. Give the dry feed either as an elemental analysis
(`elements_wt_percent`) or as relative molar amounts of registered species
(`species_moles`).

    elements_wt_percent   elemental mass percentages (ultimate analysis)
    basis                 "dry": per 100 g dry feed including ash;
                          "dry_ash_free": must sum to 100, scaled by (1 - ash/100)
    ash_wt_percent        dry-basis ash, outside the reactive model (bookkeeping only)
    unmodelled_wt_percent other dry mass outside the model
                          If either of these two is given, elements + ash + unmodelled
                          must sum to 100 (within 0.05).
    species_moles         e.g. {"cellulose" => 1.0}; replaces the elemental analysis
    moisture_wt_percent   feed moisture, % of as-received (wet) feed
    water_g_per_g_dry     added process water [g water / g dry feed]
    salts                 {name => Salt}, in mol salt / g dry feed
    species_amounts       format-2 species amounts {name => {value, unit}} including
                          water; converted with molar masses by `resolved`
"""
mutable struct Feed
    elements_wt_percent::Any
    basis::Any
    ash_wt_percent::Any
    unmodelled_wt_percent::Any
    species_moles::Any
    moisture_wt_percent::Any
    water_g_per_g_dry::Any
    salts::OrderedDict{String,Salt}
    species_amounts::Any
end
function Feed(; elements_wt_percent=nothing, basis="dry", ash_wt_percent=nothing,
              unmodelled_wt_percent=nothing, species_moles=nothing, moisture_wt_percent=0.0,
              water_g_per_g_dry=0.0, salts=OrderedDict{String,Salt}(), species_amounts=nothing)
    return Feed(elements_wt_percent, basis, ash_wt_percent, unmodelled_wt_percent, species_moles,
                moisture_wt_percent, water_g_per_g_dry, _salts(salts), species_amounts)
end
Base.:(==)(a::Feed, b::Feed) = all(getfield(a, f) == getfield(b, f) for f in fieldnames(Feed))

const FEED_FIELDS = Tuple(string.(fieldnames(Feed)))

# Feed(**data) of Python, for a dict of keyword arguments.
function _feed_from_kwargs(data::AbstractDict)
    for k in keys(data)
        string(k) in FEED_FIELDS || throw(_unexpected("Feed", k))
    end
    return Feed(; (Symbol(k) => v for (k, v) in data)...)
end

"""
    resolved(feed, molar_mass) -> Feed

This feed with species amounts converted to species_moles and a water ratio,
using `molar_mass(name) -> g/mol`. Other feeds are returned as is.
"""
function resolved(feed::Feed, molar_mass)
    feed.species_amounts === nothing && return feed
    moles, water = resolve_species(feed.species_amounts, molar_mass)
    return Feed(species_moles=moles, water_g_per_g_dry=water, salts=deepcopy(feed.salts))
end

function _unresolved(feed::Feed)
    feed.species_amounts === nothing ||
        throw(CaseError("The feed is given as species amounts; convert it with " *
                        "Feed.resolved(molar_mass) first."))
end

"""Water per dry feed [g/g]: added water plus moisture."""
function total_water_g_per_g_dry(feed::Feed)
    _unresolved(feed)
    moisture = _number(feed.moisture_wt_percent, "feed moisture_wt_percent", 0)
    moisture >= 100 && throw(CaseError("Feed moisture must be below 100 wt%."))
    return _number(feed.water_g_per_g_dry, "feed water_g_per_g_dry", 0) + moisture / (100 - moisture)
end

"""
    dry_element_wt_percent(feed, species_to_wt_percent) -> OrderedDict{String,Float64}

Dry-basis elemental wt% used by the model. `species_to_wt_percent` converts
`species_moles` (it needs the registry); `nothing` for an elemental feed.
"""
function dry_element_wt_percent(feed::Feed, species_to_wt_percent)
    _unresolved(feed)
    ash = feed.ash_wt_percent === nothing ? nothing : _number(feed.ash_wt_percent, "ash_wt_percent", 0, 100)
    other = feed.unmodelled_wt_percent === nothing ? nothing :
            _number(feed.unmodelled_wt_percent, "unmodelled_wt_percent", 0, 100)
    if feed.species_moles !== nothing
        if py_truth(ash) || py_truth(other)
            throw(CaseError("For a species-defined feed, leave ash_wt_percent and " *
                            "unmodelled_wt_percent unset or zero."))
        end
        return species_to_wt_percent(feed.species_moles)
    end
    feed.elements_wt_percent === nothing &&
        throw(CaseError("Give the feed as elements_wt_percent or species_moles."))
    analysis = OrderedDict{String,Float64}(string(e) => _number(v, "feed $e wt%", 0, 100)
                                           for (e, v) in feed.elements_wt_percent)
    if feed.basis == "dry"
        dry = analysis
    elseif feed.basis == "dry_ash_free"
        py_truth(other) && throw(CaseError("A dry_ash_free analysis requires unmodelled_wt_percent = 0."))
        py_isclose(pysum(values(analysis)), 100; abs_tol=0.05) ||
            throw(CaseError("A dry_ash_free elemental analysis must sum to 100 %."))
        a = py_truth(ash) ? ash : 0
        dry = OrderedDict{String,Float64}(e => v * (1 - a / 100) for (e, v) in analysis)
    else
        throw(CaseError("Feed basis must be 'dry' or 'dry_ash_free'."))
    end
    if ash !== nothing || other !== nothing
        total = pysum(values(dry)) + (py_truth(ash) ? ash : 0) + (py_truth(other) ? other : 0)
        if !py_isclose(total, 100; abs_tol=0.05)
            throw(CaseError("Dry elements + ash + unmodelled = $(py_format_g(total, 4)) %, not 100 %. " *
                            "No automatic renormalization is applied."))
        end
    end
    return dry
end

# ---------------------------------------------------------------- Case
# Field names and defaults, in the order of the Python dataclass. `_REQUIRED`
# marks the four fields without a default.
const _REQUIRED = :__required__
const CASE_FIELDS = (
    :temperature_C => _REQUIRED, :pressure_bar => _REQUIRED, :feed => _REQUIRED, :phases => _REQUIRED,
    :name => "case", :description => "",
    # Atomic-weight overrides [g/mol]; change feed and balance arithmetic only.
    :atomic_weights => () -> _Dict(),
    # Aqueous phase
    :aqueous_model => "electrolyte",   # "electrolyte" (ions allowed) or "molecular"
    :aqueous_enabled => true,          # false removes the aqueous phase at every temperature
    :aqueous_cutoff_C => nothing,      # aqueous phase removed at T >= cutoff (a modelling rule)
    :aqueous_water_floor => 1e-8,      # minimum aqueous water [mol/mol dry-feed C]
    # Constraints. Carbon shares count carbon atoms, relative to ALL feed carbon
    # (including salt carbon). carbon_bounds: {phase => {"exact" => f} or
    # {"min" => f, "max" => f}} (fractions 0-1); "solid" means solid carbon only
    # or all solids (solid_carbon_target).
    :carbon_bounds => () -> _Dict(),
    :solid_carbon_target => "char",    # "char" (solid carbon only) or "all_solids"
    :minimum_amounts => () -> _Dict(), # {(phase, species) => mol/mol C}
    # Old (format 1) carbon settings; moved into carbon_bounds on creation, then nothing.
    :gas_carbon_fraction => nothing, :gas_carbon_min_fraction => nothing,
    :solid_carbon_fraction => nothing, :phase_carbon_fractions => nothing,
    # Starting point: "auto", "liquid", "vapor" or "lp"
    :initialization => "auto",
    :initial_gas_water_fraction => 1e-6,  # cap on starting gas water for liquid starts
    :initial_amounts => nothing,          # {(phase, species) => mol/mol C}, replaces the LP start
    # Solver
    :solver => "ipopt",                   # "ipopt", "slsqp", "trust-constr" or "slsqp-fd"
    :starts => 1,                         # starting points per point in total
    :distinct_composition_tol => 1e-4,    # results differ above this |delta amount| [mol/mol C]
    :distinct_gibbs_tol => 1.0,           # ... or above this |delta G| [J/mol C]
    :seed => 0,
    :stability_test => nothing,           # tangent-plane test of absent phases; default: starts > 1
    :max_iterations => 5000,
    :energy_scale => nothing,             # objective scale [J/mol C]; default max(1e4, |G_start|)
    :fd_ftol => 1e-9, :fd_eps => 1e-7,    # "slsqp-fd" only
    :solver_options => () -> _Dict(),     # backend-specific, e.g. ipopt_options
    :multistart => nothing,               # old (format 1): moved into starts on creation
    # Sweeps
    :continuation => true,                # start each point from the previous solution
    :refine => nothing,                   # {start_C, stop_C, step_C}: extra temperatures in a grid
    :points => nothing,                   # explicit [[T_C, P_bar], ...], replaces the grid
    # Data to plot next to the results: {file, quantity, label}
    :comparison => nothing,
    # The feed as written in a format-2 case file; `feed` is derived from it.
    :feed_v2 => nothing,
)

"""
    Case(; temperature_C, pressure_bar, feed, phases, name="case", ...)
    Case(temperature_C, pressure_bar, feed, phases; ...)

One equilibrium calculation, or a sweep. Fields and defaults are those of the
Python dataclass (see `CASE_FIELDS`). temperature_C and pressure_bar may each
be a number or a list; lists give a full grid (every temperature at every
pressure). Alternatively `points` lists explicit [temperature_C, pressure_bar]
pairs (then set the other two to `nothing`).
"""
mutable struct Case
    temperature_C::Any
    pressure_bar::Any
    feed::Any
    phases::Any
    name::Any
    description::Any
    atomic_weights::Any
    aqueous_model::Any
    aqueous_enabled::Any
    aqueous_cutoff_C::Any
    aqueous_water_floor::Any
    carbon_bounds::Any
    solid_carbon_target::Any
    minimum_amounts::Any
    gas_carbon_fraction::Any
    gas_carbon_min_fraction::Any
    solid_carbon_fraction::Any
    phase_carbon_fractions::Any
    initialization::Any
    initial_gas_water_fraction::Any
    initial_amounts::Any
    solver::Any
    starts::Any
    distinct_composition_tol::Any
    distinct_gibbs_tol::Any
    seed::Any
    stability_test::Any
    max_iterations::Any
    energy_scale::Any
    fd_ftol::Any
    fd_eps::Any
    solver_options::Any
    multistart::Any
    continuation::Any
    refine::Any
    points::Any
    comparison::Any
    feed_v2::Any
    Case(::Val{:fields}, values::Vector{Any}) = new(values...)
end
@assert fieldnames(Case) == Tuple(first.(CASE_FIELDS))

function Case(; kwargs...)
    for k in keys(kwargs)
        hasfield(Case, k) || throw(_unexpected("Case", k))
    end
    missing_fields = [string(f) for (f, d) in CASE_FIELDS if d === _REQUIRED && !haskey(kwargs, f)]
    isempty(missing_fields) || throw(_missing_arguments("Case", missing_fields))
    values = Any[haskey(kwargs, f) ? kwargs[f] : (d isa Function ? d() : d) for (f, d) in CASE_FIELDS]
    case = Case(Val(:fields), values)
    _post_init!(case)
    return case
end
Case(temperature_C, pressure_bar, feed, phases; kwargs...) =
    Case(; temperature_C, pressure_bar, feed, phases, kwargs...)

Base.:(==)(a::Case, b::Case) = all(getfield(a, f) == getfield(b, f) for f in fieldnames(Case))

# list(x) of Python for the containers a case holds.
_pylist(x::Union{AbstractVector,Tuple}) = all(v -> v isa String, x) && !isempty(x) ? String[x...] : Any[x...]
_pylist(x::AbstractString) = Any[string(c) for c in x]
_pylist(x::AbstractDict) = Any[keys(x)...]

function _post_init!(case::Case)
    if case.feed === nothing && case.feed_v2 !== nothing
        case.feed = feed_v2_to_model(case.feed_v2)
    elseif case.feed isa AbstractDict
        case.feed = _feed_from_kwargs(case.feed)
    end
    _move_legacy_settings!(case)
    case.minimum_amounts = _pairs(case.minimum_amounts)
    case.initial_amounts === nothing || (case.initial_amounts = _pairs(case.initial_amounts))
    phases = _Dict()
    for (phase, names) in case.phases
        names isa Union{AbstractVector,Tuple,AbstractString,AbstractDict} ||
            throw(ArgumentError("$(py_repr(names)) is not iterable (phase $(py_repr(string(phase))))"))
        phases[string(phase)] = _pylist(names)
    end
    case.phases = phases
    if case.points !== nothing
        case.points = Any[_point_entry(p) for p in case.points]
    end
    return case
end

# Python's list(p) fails for an entry that is not a list; Julia gives the
# message of _explicit_points at the same moment (see README, differences).
function _point_entry(p)
    p isa Union{AbstractVector,Tuple} && return Any[p...]
    p isa AbstractString && return Any[string(c) for c in p]
    throw(CaseError("Each entry of points must be [temperature_C, pressure_bar], not $(py_repr(p))."))
end

"""Old carbon settings -> carbon_bounds; multistart -> starts (the one place
where format-1 settings are translated; the calculation reads only the new ones)."""
function _move_legacy_settings!(case::Case)
    legacy = OrderedDict{Symbol,Any}(Symbol(k) => getfield(case, Symbol(k)) for k in LEGACY_CARBON_KEYS)
    if any(v !== nothing && !(v isa AbstractDict && isempty(v)) for v in values(legacy))
        bounds, case.solid_carbon_target = legacy_carbon_bounds(; legacy..., solid_carbon_target=case.solid_carbon_target)
        case.carbon_bounds = merge_carbon_bounds(case.carbon_bounds, bounds)
    end
    for key in LEGACY_CARBON_KEYS
        setfield!(case, Symbol(key), nothing)
    end
    case.carbon_bounds = normalize_carbon_bounds(case.carbon_bounds)
    if case.multistart !== nothing
        case.starts != 1 && throw(CaseError("Give either starts or the old multistart, not both."))
        temps = try
            temperatures(case)
        catch e
            e isa Union{ValueError,ArgumentError,KeyError,MethodError} || rethrow()
            nothing
        end
        case.starts = starts_from_multistart(
            case.multistart, fluid_phase_count(case.phases, case.aqueous_enabled, case.aqueous_cutoff_C, temps))
        case.multistart = nothing
    end
    return case
end

# ---------------------------------------------------------------- grid
# sorted(set(xs)) of Python for floats: the first of equal values (0.0 == -0.0) is kept.
function _sorted_set(xs)
    out, seen = Float64[], Set{Float64}()
    for x in xs
        key = x == 0 ? 0.0 : x
        key in seen && continue
        push!(seen, key)
        push!(out, x)
    end
    return sort!(out)
end

_aslist(x) = x isa Union{AbstractVector,Tuple} ? x : Any[x]

"""
    temperatures(case) -> Vector{Float64}

Sorted temperatures [°C] of the case. In a grid this includes the refinement
points and the aqueous cut-off temperature when it lies inside the range.
"""
function temperatures(case::Case)
    case.points === nothing || return _sorted_set(t for (t, _) in _explicit_points(case))
    case.temperature_C === nothing && throw(CaseError("Give temperature_C (or points)."))
    grid = Any[_number(t, "temperature_C", -273.149999) for t in _aslist(case.temperature_C)]
    isempty(grid) && throw(CaseError("Give at least one temperature."))
    grid = refined_temperatures(grid, case.refine)
    if case.aqueous_cutoff_C !== nothing
        cutoff = py_float(case.aqueous_cutoff_C)
        minimum(grid) < cutoff < maximum(grid) && push!(grid, cutoff)
    end
    return _sorted_set(py_round(py_float(t), 8) for t in grid)
end

"""Sorted pressures [bar] of the case."""
function pressures(case::Case)
    case.points === nothing || return _sorted_set(p for (_, p) in _explicit_points(case))
    case.pressure_bar === nothing && throw(CaseError("Give pressure_bar (or points)."))
    values = _aslist(case.pressure_bar)
    isempty(values) && throw(CaseError("Give at least one pressure."))
    return _sorted_set(py_round(_number(p, "pressure_bar", 1e-300), 8) for p in values)
end

function _explicit_points(case::Case)
    out = Tuple{Float64,Float64}[]
    for item in case.points
        if !(item isa Union{AbstractVector,Tuple}) || length(item) != 2
            throw(CaseError("Each entry of points must be [temperature_C, pressure_bar], not $(py_repr(item))."))
        end
        point = (py_round(_number(item[1], "points temperature_C", -273.149999), 8),
                 py_round(_number(item[2], "points pressure_bar", 1e-300), 8))
        point in out || push!(out, point)
    end
    isempty(out) && throw(CaseError("points is empty."))
    return out
end

"""
    grid(case) -> Vector{Tuple{Float64,Float64}}

All (temperature_C, pressure_bar) points: explicit points in their order, else
every temperature at every pressure (by pressure, then temperature).
"""
function grid(case::Case)
    case.points === nothing || return _explicit_points(case)
    temps = temperatures(case)
    return [(T, P) for P in pressures(case) for T in temps]
end

"""
    chains(case) -> Vector{Vector{Tuple{Float64,Float64}}}

Groups of points that must be solved in order (continuation). Without
continuation every point is its own group. With continuation a grid gives one
chain per pressure (rising temperature) and explicit points one chain in list order.
"""
function chains(case::Case)
    points = grid(case)
    py_truth(case.continuation) || return [[point] for point in points]
    case.points === nothing || return [points]
    out = OrderedDict{Float64,Vector{Tuple{Float64,Float64}}}()
    for (T, P) in points
        push!(get!(out, P, Tuple{Float64,Float64}[]), (T, P))
    end
    return collect(values(out))
end

is_single_point(case::Case) = length(grid(case)) == 1

"""
    at(case, temperature_C, pressure_bar=nothing) -> Case

Copy of the case at one temperature (and pressure; optional when the case has
a single pressure).
"""
function at(case::Case, temperature_C, pressure_bar=nothing)
    if pressure_bar === nothing
        ps = pressures(case)
        length(ps) == 1 || throw(CaseError("This case has several pressures; give pressure_bar too."))
        pressure_bar = ps[1]
    end
    point = deepcopy(case)
    point.temperature_C = py_float(temperature_C)
    point.pressure_bar = py_float(pressure_bar)
    point.points = nothing
    point.refine = nothing
    return point
end

"""
    replace(case; changes...) -> Case

Copy of the case with some fields changed, normalised again (Python:
`case.replace(**changes)`). A new `feed` drops `feed_v2`; a new `feed_v2`
gives a new `feed`.
"""
function Base.replace(case::Case; changes...)
    point = deepcopy(case)
    for (key, value) in changes
        hasfield(Case, key) || throw(ArgumentError("Case has no field $(py_repr(string(key)))"))
        setfield!(point, key, value)
    end
    if haskey(changes, :feed) && !haskey(changes, :feed_v2)
        point.feed_v2 = nothing  # the written feed no longer describes the new one
    elseif haskey(changes, :feed_v2) && !haskey(changes, :feed) && point.feed_v2 !== nothing
        point.feed = feed_v2_to_model(point.feed_v2; salts=case.feed === nothing ? nothing : case.feed.salts)
    end
    return _post_init!(point)
end

"""
    validate(case) -> case

Checks that need no databases. Throws a `CaseError` with a readable message.
Like Python it also stores `starts` as an Int.
"""
function validate(case::Case)
    unknown = setdiff(keys(case.phases), PHASES)
    isempty(unknown) || throw(CaseError("Unknown phase(s) $(py_repr(py_sorted(unknown))); use $(py_repr(collect(PHASES)))."))
    case.solver in SOLVERS || throw(CaseError("solver must be one of $(py_repr(SOLVERS)), not $(py_repr(case.solver))."))
    case.initialization in INITIALIZATIONS || throw(CaseError("initialization must be one of $(py_repr(INITIALIZATIONS))."))
    case.aqueous_model in ("electrolyte", "molecular") ||
        throw(CaseError("aqueous_model must be 'electrolyte' or 'molecular'."))
    case.solid_carbon_target in ("char", "all_solids") ||
        throw(CaseError("solid_carbon_target must be 'char' or 'all_solids'."))
    if case.points !== nothing && py_truth(case.refine)
        throw(CaseError("refine applies to a temperature grid; it cannot be used with points."))
    end
    check_carbon_bounds(case.carbon_bounds)
    s = case.starts
    if s isa Real && !(s isa Bool) && !isfinite(s)
        # Python: int(nan) / int(inf) fail inside the check with these messages.
        throw(CaseError(isnan(s) ? "cannot convert float NaN to integer" : "cannot convert float infinity to integer"))
    end
    if s isa Bool || !(s isa Real) || s != trunc(s) || s < 1
        throw(CaseError("starts (initial guesses) must be a whole number of at least 1, not $(py_repr(s))."))
    end
    case.starts = Int(s)
    for name in (:distinct_composition_tol, :distinct_gibbs_tol)
        _number(getfield(case, name), string(name), 0)
    end
    grid(case)
    return case
end

function _number(value, name, low=nothing, high=nothing)
    v = try
        py_float(value)
    catch e
        e isa ArgumentError || rethrow()
        throw(CaseError("$name must be a number, not $(py_repr(value))."))
    end
    if !isfinite(v) || (low !== nothing && v < low) || (high !== nothing && v > high)
        throw(CaseError("Invalid $name: $(py_str(v))"))
    end
    return v
end

"""{(phase, species) => v} from either that form or nested {phase => {species => v}}."""
function _pairs(mapping)
    out = OrderedDict{Tuple{String,String},Float64}()
    mapping === nothing && return out
    for (key, value) in mapping
        if key isa Tuple
            out[(string(key[1]), string(key[2]))] = py_float(value)
        else
            for (species, amount) in value
                out[(string(key), string(species))] = py_float(amount)
            end
        end
    end
    return out
end

function _nested(pairs)
    out = _Dict()
    for ((phase, species), value) in pairs
        get!(out, phase, _Dict())[species] = value
    end
    return out
end

# ---------------------------------------------------------------- carbon bounds
const PHASE_LABELS = OrderedDict("gas" => "Gas", "aqueous" => "Aqueous", "organic_liquid" => "Organic liquid",
                                 "solid" => "Solid")
const FLUID_PHASES = ("gas", "aqueous", "organic_liquid")
const BOUND_KINDS = ("exact", "min", "max")
const LEGACY_CARBON_KEYS = ("gas_carbon_fraction", "gas_carbon_min_fraction", "solid_carbon_fraction",
                            "phase_carbon_fractions")

_label(phase) = get(PHASE_LABELS, phase, py_str(phase))

"""
    legacy_carbon_bounds(; gas_carbon_fraction, gas_carbon_min_fraction,
                         solid_carbon_fraction, phase_carbon_fractions,
                         solid_carbon_target="char") -> (carbon_bounds, solid_carbon_target)

gas_carbon_fraction -> gas.exact, gas_carbon_min_fraction -> gas.min,
solid_carbon_fraction -> solid.exact, phase_carbon_fractions[p] -> p.exact.
The old phase_carbon_fractions["solid"] counted all solids, so it sets
solid_carbon_target = "all_solids" (it cannot be combined with
solid_carbon_fraction, which counted char only by default).
"""
function legacy_carbon_bounds(; gas_carbon_fraction=nothing, gas_carbon_min_fraction=nothing,
                              solid_carbon_fraction=nothing, phase_carbon_fractions=nothing,
                              solid_carbon_target="char")
    bounds = _Dict()
    if gas_carbon_fraction !== nothing && gas_carbon_min_fraction !== nothing
        throw(CaseError("Use either gas_carbon_fraction (equality) or " *
                        "gas_carbon_min_fraction (lower bound), not both"))
    end
    gas_carbon_fraction === nothing || (bounds["gas"] = _Dict("exact" => gas_carbon_fraction))
    gas_carbon_min_fraction === nothing || (bounds["gas"] = _Dict("min" => gas_carbon_min_fraction))
    solid_carbon_fraction === nothing || (bounds["solid"] = _Dict("exact" => solid_carbon_fraction))
    for (phase, fraction) in (py_truth(phase_carbon_fractions) ? phase_carbon_fractions : ())
        if phase == "gas" && gas_carbon_fraction !== nothing
            if !py_isclose(py_float(fraction), py_float(gas_carbon_fraction); rel_tol=0, abs_tol=1e-10)
                throw(CaseError("Conflicting gas-carbon constraints"))
            end
            continue
        end
        phase == "gas" && gas_carbon_min_fraction !== nothing && throw(CaseError("Conflicting gas-carbon constraints"))
        if phase == "solid"
            if solid_carbon_fraction !== nothing
                throw(CaseError("phase_carbon_fractions.solid (all solids) cannot be combined with " *
                                "solid_carbon_fraction; use [constraints.carbon] solid and " *
                                "solid_carbon_target."))
            end
            solid_carbon_target = "all_solids"
        end
        bounds[phase] = _Dict("exact" => fraction)
    end
    return bounds, solid_carbon_target
end

function merge_carbon_bounds(bounds, extra)
    out = _Dict()
    for (phase, b) in (py_truth(bounds) ? bounds : ())
        out[phase] = b isa AbstractDict ? _Dict(b) : b
    end
    for (phase, bound) in extra
        if haskey(out, phase) && out[phase] != bound
            throw(CaseError("$(_label(phase)): the carbon share is given twice (old and new " *
                            "keys) with different values."))
        end
        out[phase] = _Dict(bound)
    end
    return out
end

"""Ordered copy: gas, solid, then the other phases as given; numbers as floats."""
function normalize_carbon_bounds(bounds)
    bounds = py_truth(bounds) ? bounds : _Dict()
    order = [[p for p in ("gas", "solid") if haskey(bounds, p)]; [p for p in keys(bounds) if !(p in ("gas", "solid"))]]
    out = _Dict()
    for phase in order
        bound = bounds[phase]
        if !(bound isa AbstractDict)
            out[phase] = bound
            continue
        end
        item = _Dict()
        for key in [[k for k in BOUND_KINDS if haskey(bound, k)]; [k for k in keys(bound) if !(k in BOUND_KINDS)]]
            value = bound[key]
            item[key] = value isa Real && !(value isa Bool) ? Float64(value) : value
        end
        out[phase] = item
    end
    return out
end

"""Throw a CaseError for the first malformed bound (needs no species data)."""
function check_carbon_bounds(bounds)
    for (phase, bound) in (py_truth(bounds) ? bounds : ())
        phase in PHASES ||
            throw(CaseError("Unknown phase $(py_repr(phase)) in [constraints.carbon]; use $(join(PHASES, ", "))."))
        label = _label(phase)
        if !(bound isa AbstractDict) || isempty(bound) || !all(k -> k in BOUND_KINDS, keys(bound))
            throw(CaseError("$label carbon: use exact, or min and/or max (e.g. $phase = { min = 0.9 })."))
        end
        if haskey(bound, "exact") && (haskey(bound, "min") || haskey(bound, "max"))
            throw(CaseError("$label carbon: exact cannot be combined with min or max."))
        end
        for (kind, value) in bound
            if value isa Bool || !(value isa Real) || !isfinite(value) || !(0 <= value <= 1)
                throw(CaseError("$label carbon: $kind must be a fraction between 0 and 1 " *
                                "(0 % to 100 %), not $(py_repr(value))."))
            end
        end
        if haskey(bound, "min") && haskey(bound, "max") && bound["min"] > bound["max"]
            throw(CaseError("$label: the minimum is larger than the maximum."))
        end
    end
end

"""
    carbon_bound_messages(case, carbon_solids=nothing) -> Vector{String}

Readable problems with the carbon bounds of `case` (empty: none found).
`carbon_solids`: optional {solid species name => true if it contains carbon};
used to decide whether "char only" covers all solid carbon.
"""
function carbon_bound_messages(case::Case, carbon_solids=nothing)
    messages = String[]
    bounds = case.carbon_bounds
    try
        check_carbon_bounds(bounds)
    catch e
        e isa CaseError || rethrow()
        return [e.msg]
    end
    present = [p for p in PHASES if py_truth(get(case.phases, p, nothing))]   # PHASES order (Python: a set)
    case.aqueous_enabled === false && filter!(!=("aqueous"), present)
    for phase in keys(bounds)
        phase in present || push!(messages, "$(_label(phase)) has a carbon bound but no species.")
    end
    fixed = pysum(get(b, "exact", 0.0) + get(b, "min", 0.0) for b in values(bounds))
    if fixed > 1 + 1e-10
        push!(messages, "The carbon shares add up to more than 100 %.")
    elseif !isempty(present) && all(p -> haskey(get(bounds, p, _Dict()), "exact"), present)
        solid_complete = !("solid" in present) || case.solid_carbon_target == "all_solids" ||
            (carbon_solids !== nothing &&
             !any(py_truth(get(carbon_solids, s, true)) for s in case.phases["solid"] if s != "carbon"))
        total = pysum(bounds[p]["exact"] for p in present)
        if solid_complete && !py_isclose(total, 1.0; rel_tol=0, abs_tol=1e-10)
            push!(messages, "With every phase fixed, the carbon shares must add up to 100 %.")
        end
    end
    return messages
end

"""
    effective_carbon_bounds(bounds, aqueous_disabled) -> (bounds, dropped)

Bounds at a point: without an aqueous phase the aqueous bound is dropped;
dropped = ["aqueous carbon"] when that bound asked for carbon there.
"""
function effective_carbon_bounds(bounds, aqueous_disabled)
    (!py_truth(aqueous_disabled) || !haskey(bounds, "aqueous")) && return _Dict(bounds), String[]
    bound = bounds["aqueous"]
    out = _Dict(p => b for (p, b) in bounds if p != "aqueous")
    asked = get(bound, "exact", 0.0) > 0 || get(bound, "min", 0.0) > 0
    return out, asked ? ["aqueous carbon"] : String[]
end

# ---------------------------------------------------------------- starting points
"""
Fluid phases (gas, aqueous, organic liquid) with at least one species. The
aqueous phase does not count when it is disabled, or removed by the cut-off at
every temperature of the case.
"""
function fluid_phase_count(phases, aqueous_enabled=true, aqueous_cutoff_C=nothing, temperatures=nothing)
    count = 0
    for phase in FLUID_PHASES
        py_truth(phases === nothing ? nothing : get(phases, phase, nothing)) || continue
        if phase == "aqueous"
            aqueous_enabled === false && continue
            if aqueous_cutoff_C !== nothing && py_truth(temperatures) &&
               minimum(temperatures) >= py_float(aqueous_cutoff_C)
                continue
            end
        end
        count += 1
    end
    return count
end

"""
Total starting points of the old multistart = k: 1 + 2 F + k (k = 0: 1). The
old code made the initial start, a phase-rich and a phase-lean LP start for
each of the F fluid phases, and k random starts.
"""
function starts_from_multistart(multistart, fluid_phases)
    k = py_int(py_truth(multistart) ? multistart : 0)
    k < 0 && throw(CaseError("multistart cannot be negative."))
    return k > 0 ? 1 + 2 * py_int(fluid_phases) + k : 1
end

# numpy.arange(start, stop, step) for floats, element by element as NumPy fills it.
function _np_arange(start::Float64, stop::Float64, step::Float64)
    n = ceil((stop - start) / step)
    n <= 0 && return Float64[]
    n = Int(n)
    out = Vector{Float64}(undef, n)
    out[1] = start
    n >= 2 && (out[2] = start + step)
    delta = n >= 2 ? out[2] - out[1] : 0.0
    for i in 3:n
        out[i] = start + (i - 1) * delta
    end
    return out
end

"""`grid` plus the extra temperatures of the old `refine` setting
({start_C, stop_C, step_C}, limited to the range of the grid)."""
function refined_temperatures(grid, refine)
    grid = Any[grid...]
    if py_truth(refine)
        step = _number(refine["step_C"], "refine step_C", 1e-12)
        a = max(py_float(minimum(grid)), py_float(refine["start_C"]))
        b = min(py_float(maximum(grid)), py_float(refine["stop_C"]))
        b >= a && append!(grid, _np_arange(a, b + step * 1e-8, step))
    end
    return grid
end

# ---------------------------------------------------------------- old (format 1) files
const NEW_KEYS = OrderedDict("feed" => ("total_kg", "solids", "ash", "species"), "constraints" => ("carbon",),
                             "solver" => ("starts", "distinct_composition_tol", "distinct_gibbs_tol"))

"""True for case-file data in format 2 (format = 2, or any key only format 2 has)."""
function is_v2(data)
    get(data, "format", nothing) == 2 && return true
    return any(get(data, section, nothing) isa AbstractDict && haskey(data[section], key)
               for (section, keys) in NEW_KEYS for key in keys)
end

"""
    legacy_to_v2(data, molar_mass=nothing) -> OrderedDict

Case-file data in format 2 for old (format-1) data, giving the same calculation.
Throws LegacyConversionError when the file cannot be represented (its message is
meant for the user): unmodelled_wt_percent > 0, or a species-moles feed without
`molar_mass(name) -> g/mol`. Data already in format 2 is returned as a copy.
See `legacy_feed_to_v2` for the feed rules.
"""
function legacy_to_v2(data, molar_mass=nothing)
    data = deepcopy(data)
    is_v2(data) && return data
    out = _Dict("format" => 2)
    for (k, v) in data
        out[k] = v
    end
    haskey(out, "feed") && (out["feed"] = legacy_feed_to_v2(out["feed"], molar_mass))
    constraints = get(out, "constraints", nothing)
    if constraints isa AbstractDict
        legacy = OrderedDict{Symbol,Any}(Symbol(key) => pop!(constraints, key) for key in LEGACY_CARBON_KEYS
                                         if haskey(constraints, key))
        if !isempty(legacy)
            target0 = get(constraints, "solid_carbon_target", "char")
            bounds, target = legacy_carbon_bounds(; legacy..., solid_carbon_target=target0)
            target != target0 && (constraints["solid_carbon_target"] = target)
            isempty(bounds) || (constraints["carbon"] = normalize_carbon_bounds(bounds))
        end
    end
    solver = get(out, "solver", nothing)
    if solver isa AbstractDict && haskey(solver, "multistart")
        conditions = get(out, "conditions", _Dict())
        aqueous = get(out, "aqueous", _Dict())
        temps = nothing
        if haskey(conditions, "temperature_C")
            values = conditions["temperature_C"]
            sweep = get(out, "sweep", nothing)
            temps = refined_temperatures(values isa AbstractVector ? values : Any[values],
                                         get(py_truth(sweep) ? sweep : _Dict(), "refine", nothing))
        elseif haskey(conditions, "points")
            temps = [p[1] for p in conditions["points"]]
        end
        count = fluid_phase_count(get(out, "phases", _Dict()), get(aqueous, "enabled", true),
                                  get(aqueous, "cutoff_C", nothing), temps)
        multistart = pop!(solver, "multistart")
        solver["starts"] = starts_from_multistart(multistart, count)
    end
    sweep = get(out, "sweep", nothing)
    if sweep isa AbstractDict && haskey(sweep, "refine")
        refine = pop!(sweep, "refine")
        conditions = get(out, "conditions", _Dict())
        if py_truth(refine) && haskey(conditions, "temperature_C") && !haskey(conditions, "points")
            values = conditions["temperature_C"]
            merged = _sorted_set(py_round(py_float(t), 8)
                                 for t in refined_temperatures(values isa AbstractVector ? values : Any[values], refine))
            conditions["temperature_C"] = Any[isinteger(t) ? Int(t) : t for t in merged]
        end
        isempty(sweep) && delete!(out, "sweep")
    end
    return out
end

# ---------------------------------------------------------------- TOML case files
# section -> {key in file => Case field}. Old (format 1) keys are still read.
const SECTIONS = OrderedDict(
    "conditions" => OrderedDict("temperature_C" => "temperature_C", "pressure_bar" => "pressure_bar",
                                "points" => "points"),
    "aqueous" => OrderedDict("model" => "aqueous_model", "enabled" => "aqueous_enabled",
                             "cutoff_C" => "aqueous_cutoff_C", "water_floor" => "aqueous_water_floor"),
    "constraints" => OrderedDict("carbon" => "carbon_bounds", "solid_carbon_target" => "solid_carbon_target",
                                 "minimum_amounts" => "minimum_amounts",
                                 (k => k for k in LEGACY_CARBON_KEYS)...),
    "initial_guess" => OrderedDict("initialization" => "initialization",
                                   "initial_gas_water_fraction" => "initial_gas_water_fraction",
                                   "amounts" => "initial_amounts"),
    "solver" => OrderedDict("method" => "solver", "starts" => "starts",
                            "distinct_composition_tol" => "distinct_composition_tol",
                            "distinct_gibbs_tol" => "distinct_gibbs_tol", "seed" => "seed",
                            "stability_test" => "stability_test", "max_iterations" => "max_iterations",
                            "energy_scale" => "energy_scale", "fd_ftol" => "fd_ftol", "fd_eps" => "fd_eps",
                            "options" => "solver_options", "multistart" => "multistart"),
    "sweep" => OrderedDict("continuation" => "continuation", "refine" => "refine"),
    "comparison" => OrderedDict("file" => "file", "quantity" => "quantity", "label" => "label"),
)
const TOP_LEVEL = ("name", "description", "atomic_weights")
const FORMAT = 2
const LEGACY_FEED_ONLY = ("basis", "ash_wt_percent", "unmodelled_wt_percent", "species_moles",
                          "moisture_wt_percent", "water_g_per_g_dry")

"""(Feed, feed_v2) for a [feed] table of either format."""
function _feed_from_dict(feed, version)
    feed isa AbstractDict || throw(CaseError("[feed] must be a table."))
    if is_v2_feed(feed)
        mixed = [k for k in LEGACY_FEED_ONLY if haskey(feed, k)]
        if !isempty(mixed)
            throw(CaseError("[feed] mixes the new keys (total_kg, solids, ash, species) with old " *
                            "ones $(py_repr(mixed)); remove the old ones."))
        end
        written = _Dict(k => deepcopy(v) for (k, v) in feed if k != "salts")
        return feed_v2_to_model(feed), written
    end
    if version == FORMAT && !any(haskey(feed, k) for k in LEGACY_FEED_ONLY)
        throw(CaseError("[feed] needs total_kg and solids (and elements_wt_percent), " *
                        "or a [feed.species] table."))
    end
    known = [f for f in FEED_FIELDS if f != "species_amounts"]
    unknown = setdiff(keys(feed), known)
    isempty(unknown) ||
        throw(CaseError("Unknown [feed] key(s) $(py_repr(py_sorted(unknown))); allowed: $(py_repr(py_sorted(known)))."))
    return _feed_from_kwargs(feed), nothing
end

"""
    case_from_dict(data; name=nothing) -> Case

Case from case-file data (nested dicts as `parse_toml` returns them), format 2
or the old format (or a mixture). Validated.
"""
function case_from_dict(data; name=nothing)
    data = _Dict(deepcopy(data))
    version = pop!(data, "format", nothing)
    if !(version === nothing || version == 1 || version == FORMAT)
        throw(CaseError("format = $(py_repr(version)) is not known; this version reads format 2 " *
                        "(and old files without format)."))
    end
    kwargs = OrderedDict{Symbol,Any}()
    for key in TOP_LEVEL
        haskey(data, key) && (kwargs[Symbol(key)] = pop!(data, key))
    end
    haskey(data, "feed") || throw(CaseError("The case needs a [feed] section."))
    kwargs[:feed], kwargs[:feed_v2] = _feed_from_dict(pop!(data, "feed"), version)
    kwargs[:phases] = pop!(data, "phases", _Dict())
    for (section, keys_) in SECTIONS
        values = pop!(data, section, _Dict())
        values isa AbstractDict || throw(CaseError("[$section] must be a table."))
        if section == "comparison"
            if !isempty(values)
                unknown = setdiff(keys(values), keys(keys_))
                if !isempty(unknown) || !haskey(values, "file")
                    throw(CaseError("[comparison] needs file (and optionally quantity, label)."))
                end
                kwargs[:comparison] = values
            end
            continue
        end
        unknown = setdiff(keys(values), keys(keys_))
        isempty(unknown) ||
            throw(CaseError("Unknown [$section] key(s) $(py_repr(py_sorted(unknown))); allowed: $(py_repr(py_sorted(keys(keys_))))."))
        for (key, value) in values
            kwargs[Symbol(keys_[key])] = value
        end
    end
    isempty(data) || throw(CaseError("Unknown section(s) or key(s): $(py_repr(py_sorted(keys(data))))."))
    if haskey(kwargs, :points)
        if haskey(kwargs, :temperature_C) || haskey(kwargs, :pressure_bar)
            throw(CaseError("[conditions] has points; leave out temperature_C and pressure_bar."))
        end
        kwargs[:temperature_C] = kwargs[:pressure_bar] = nothing
    else
        for required in (:temperature_C, :pressure_bar)
            haskey(kwargs, required) || throw(CaseError("[conditions] needs $required (or points)."))
        end
    end
    if py_truth(name) && !haskey(kwargs, :name)
        kwargs[:name] = name
    end
    return validate(Case(; kwargs...))
end

"""
    load_case(path) -> Case

Read a TOML case file (key order kept). The comparison file becomes an
absolute path, relative to the case file.
"""
function load_case(path::AbstractString)
    stem = splitext(basename(path))[1]
    case = case_from_dict(read_toml(path); name=stem)
    if py_truth(case.comparison)
        file = abspath(joinpath(dirname(abspath(path)), case.comparison["file"]))
        case.comparison["file"] = ispath(file) ? realpath(file) : file
    end
    return case
end
