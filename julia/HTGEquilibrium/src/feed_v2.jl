# The feed of case-file format 2, as plain functions (Python: htg/feed_v2.py).
#
# Elemental analysis:
#
#     [feed]
#     total_kg = 10.0                          # solids + water
#     solids = { value = 15.0, unit = "%" }    # "%" of the total, or "kg"
#     ash = { value = 8.45, unit = "%" }       # "%" of the solids, or "kg"; omit = no ash
#     elements_wt_percent = { C = 30.4, H = 4.41, O = 53.58, N = 2.66, S = 0.51 }
#
# The element percentages are of the ash-free solids, as entered; they are
# scaled to sum to 100 when used.
#
# Species amounts:
#
#     [feed.species]
#     water = { value = 4.5, unit = "kg" }     # always present; "kg" or "mol"
#     cellulose = { value = 0.5, unit = "kg" }
#
# `feed_v2_to_model` maps either form onto the model's `Feed`; `feed_summary`
# gives the numbers the Case tab shows. Species need molar masses: pass
# `molar_mass(name) -> g/mol` (see `molar_mass_function`) or
# `species_elements(name) -> {element => count}`. Python's
# `registry_elements(registry)` becomes `species_elements_function(entries)`
# with the registry's {name => elements} table (the registry is ported separately).

"""The feed cannot be used; the message names the field (Python: FeedError, a ValueError)."""
struct FeedError <: ValueError
    msg::String
end

const ELEMENTAL_KEYS = ("total_kg", "solids", "ash", "elements_wt_percent")
const SPECIES_KEY = "species"
const V2_FEED_KEYS = ("total_kg", "solids", "ash", "species")  # keys that only exist in format 2
const SOLIDS_UNITS = ("%", "kg")
const SPECIES_UNITS = ("kg", "mol")

# ---------------------------------------------------------------- helpers
function _feed_number(value, field)
    value isa Bool && throw(FeedError("$field must be a number, not $(py_repr(value))."))
    number = try
        py_float(value)
    catch e
        e isa ArgumentError || rethrow()
        throw(FeedError("$field must be a number, not $(py_repr(value))."))
    end
    isfinite(number) || throw(FeedError("$field must be a finite number."))
    return number
end

# a / b, with Python's ZeroDivisionError (here DivideError) for b == 0.
_pydiv(a, b) = b == 0 ? throw(DivideError()) : a / b

"""(value, unit) of an entry `{value = x, unit = "u"}`."""
function amount(entry, field, units)
    if !(entry isa AbstractDict) || !haskey(entry, "value")
        throw(FeedError("$field: write { value = ..., unit = \"$(units[1])\" }."))
    end
    unknown = setdiff(keys(entry), ("value", "unit"))
    isempty(unknown) ||
        throw(FeedError("$field: unknown key(s) $(py_repr(py_sorted(unknown))); use value and unit."))
    unit = get(entry, "unit", units[1])
    if !(unit in units)
        throw(FeedError("$field: the unit must be " * join(("\"$u\"" for u in units), " or ") *
                        ", not $(py_repr(unit))."))
    end
    return _feed_number(entry["value"], field), unit
end

"""\"species\" or \"elements\" for a format-2 [feed] table."""
feed_mode(feed) = haskey(feed, SPECIES_KEY) ? "species" : "elements"

is_v2_feed(feed) = feed isa AbstractDict && any(haskey(feed, key) for key in V2_FEED_KEYS)

function _atomic_weights(atomic_weights)
    aw = Dict{String,Any}(ATOMIC_WEIGHTS)
    py_truth(atomic_weights) && for (e, w) in atomic_weights
        aw[string(e)] = w
    end
    return aw
end

water_molar_mass(atomic_weights=nothing) = (aw = _atomic_weights(atomic_weights); 2.0 * aw["H"] + aw["O"])

"""
    molar_mass_function(species_elements, atomic_weights=nothing) -> molar_mass

`molar_mass(name)` [g/mol] from `species_elements(name) -> {element => count}`.
"""
function molar_mass_function(species_elements, atomic_weights=nothing)
    aw = _atomic_weights(atomic_weights)
    function molar_mass(name)
        elements = species_elements(name)
        missing_ = [e for e in keys(elements) if !haskey(aw, e)]
        isempty(missing_) ||
            throw(FeedError("No atomic weight for $(join(missing_, ", ")) (species $(py_str(name)))."))
        return Float64(sum(aw[e] * count for (e, count) in elements; init=0.0))
    end
    return molar_mass
end

"""
    species_elements_function(entries) -> species_elements

`species_elements(name) -> {element => count}` from a table
{name => elements} (or {name => entry with "elements"}). For a `Registry`
use `registry_elements`.
"""
function species_elements_function(entries::AbstractDict)
    function species_elements(name)
        haskey(entries, name) || throw(FeedError("Unknown feed species $(py_repr(name))."))
        entry = entries[name]
        elements = entry isa AbstractDict && haskey(entry, "elements") ? entry["elements"] : entry
        return OrderedDict{String,Any}(elements)
    end
    return species_elements
end

"""
    registry_elements(registry) -> species_elements

`species_elements(name)` for a species `Registry` (molecules, solids and ions;
later tables win, as Python's `{**all_molecular, **solid, **ion}`).
"""
registry_elements(registry) = species_elements_function(merge(all_molecular(registry), registry.solid, registry.ion))

# ---------------------------------------------------------------- elemental analysis
"""(used, entered sum): the element percentages scaled to sum to 100."""
function scaled_elements(elements)
    (!(elements isa AbstractDict) || isempty(elements)) && throw(FeedError("Elements: give at least one element."))
    entered = OrderedDict{String,Float64}()
    for (element, value) in elements
        number = _feed_number(value, "Element $element")
        number < 0 && throw(FeedError("Element $element: the percentage cannot be negative."))
        entered[string(element)] = number
    end
    total = sum(values(entered); init=0.0)
    total <= 0 && throw(FeedError("Elements: the percentages add up to zero."))
    abs(total - 100.0) <= 1e-12 * 100.0 && return copy(entered), total  # already 100: keep exactly
    return OrderedDict{String,Float64}(e => v * 100.0 / total for (e, v) in entered), total
end

"""Masses [kg] of an elemental-analysis feed: total, solids, water, ash, ash_free."""
function elemental_masses(feed)
    haskey(feed, "total_kg") || throw(FeedError("Total mass: give total_kg (solids + water)."))
    total = _feed_number(feed["total_kg"], "Total mass")
    total <= 0 && throw(FeedError("Total mass must be positive."))
    haskey(feed, "solids") || throw(FeedError("Solids: give solids = { value = ..., unit = \"%\" or \"kg\" }."))
    value, unit = amount(feed["solids"], "Solids", SOLIDS_UNITS)
    solids = unit == "%" ? value / 100.0 * total : value
    solids <= 0 && throw(FeedError("Solids must be positive."))
    solids > total * (1 + 1e-12) && throw(FeedError("Solids can't be more than the total mass."))
    solids = total < solids ? total : solids     # Python min(solids, total)
    ash = 0.0
    if get(feed, "ash", nothing) !== nothing
        value, unit = amount(feed["ash"], "Ash", SOLIDS_UNITS)
        value < 0 && throw(FeedError("Ash cannot be negative."))
        ash = unit == "%" ? value / 100.0 * solids : value
        ash > solids * (1 + 1e-12) && throw(FeedError("Ash can't be more than the solids."))
    end
    return OrderedDict{String,Float64}("total" => total, "solids" => solids, "water" => total - solids,
                                       "ash" => ash, "ash_free" => solids - ash)
end

"""Ash as % of the solids, or `nothing` when the feed gives no ash."""
function ash_percent(feed)
    get(feed, "ash", nothing) === nothing && return nothing
    value, unit = amount(feed["ash"], "Ash", SOLIDS_UNITS)
    unit == "%" && return value
    return 100.0 * elemental_masses(feed)["ash"] / elemental_masses(feed)["solids"]
end

# ---------------------------------------------------------------- species amounts
"""{name => (kg, mol)} of a species-amounts feed (water included)."""
function species_masses(species, molar_mass)
    if !(species isa AbstractDict) || !haskey(species, "water")
        throw(FeedError("Species amounts: water is always present; give water = { value = ..., unit = \"kg\" }."))
    end
    out = OrderedDict{String,Tuple{Float64,Float64}}()
    for (name, entry) in species
        value, unit = amount(entry, "Species $name", SPECIES_UNITS)
        value < 0 && throw(FeedError("Species $name: the amount cannot be negative."))
        M = molar_mass(name)
        out[name] = unit == "kg" ? (value, _pydiv(1000.0 * value, M)) : (value * M / 1000.0, value)
    end
    if !any(kg > 0 for (name, (kg, _)) in out if name != "water")
        throw(FeedError("Species amounts: add at least one species besides water."))
    end
    return out
end

"""(species_moles, water_g_per_g_dry) of a species-amounts feed."""
function resolve_species(species, molar_mass)
    masses = species_masses(species, molar_mass)
    moles = OrderedDict{String,Float64}(name => mol for (name, (_, mol)) in masses if name != "water")
    dry_kg = sum(kg for (name, (kg, _)) in masses if name != "water"; init=0.0)
    return moles, masses["water"][1] / dry_kg
end

# ---------------------------------------------------------------- mapping to the model
"""
    feed_v2_to_model(feed; salts=nothing, molar_mass=nothing) -> Feed

The model `Feed` for a format-2 [feed] table (`salts`: {name => Salt or dict}).
Elemental analysis: elements scaled to 100 on the dry, ash-free basis,
ash_wt_percent = 100 ash / solids, water_g_per_g_dry = water / solids.
Species amounts: species_moles and water_g_per_g_dry = water mass / non-water
mass; without `molar_mass` the conversion is left to the calculation
(`Feed.species_amounts`, see `resolved`). Positional `salts, molar_mass` work too.
"""
function feed_v2_to_model(feed, salts_=nothing, molar_mass_=nothing; salts=salts_, molar_mass=molar_mass_)
    salt_table = py_truth(salts) ? salts : (py_truth(get(feed, "salts", nothing)) ? feed["salts"] : _Dict())
    if feed_mode(feed) == "species"
        extra = setdiff(keys(feed), (SPECIES_KEY, "salts"))
        isempty(extra) || throw(FeedError("[feed.species] replaces $(py_repr(py_sorted(extra))); remove them."))
        if molar_mass === nothing
            species_masses(feed[SPECIES_KEY], name -> 1.0)  # checks the entries
            return Feed(species_amounts=deepcopy(feed[SPECIES_KEY]), salts=salt_table)
        end
        moles, water = resolve_species(feed[SPECIES_KEY], molar_mass)
        return Feed(species_moles=moles, water_g_per_g_dry=water, salts=salt_table)
    end
    unknown = setdiff(keys(feed), ELEMENTAL_KEYS, ("salts",))
    if !isempty(unknown)
        throw(FeedError("Unknown [feed] key(s) $(py_repr(py_sorted(unknown))); allowed: " *
                        "$(py_repr(collect(ELEMENTAL_KEYS))) (or [feed.species]) and [feed.salts]."))
    end
    masses = elemental_masses(feed)
    haskey(feed, "elements_wt_percent") ||
        throw(FeedError("Elements: give elements_wt_percent (mass % of the ash-free solids)."))
    used, _ = scaled_elements(feed["elements_wt_percent"])
    ash = get(feed, "ash", nothing) === nothing ? nothing : 100.0 * masses["ash"] / masses["solids"]
    if ash !== nothing && get(feed["ash"], "unit", "%") == "%"
        ash = _feed_number(feed["ash"]["value"], "Ash")  # exactly as written
    end
    return Feed(elements_wt_percent=used, basis="dry_ash_free", ash_wt_percent=ash, moisture_wt_percent=0.0,
                water_g_per_g_dry=masses["water"] / masses["solids"], salts=salt_table)
end

# ---------------------------------------------------------------- summary for the form
_row(name, kg, mol) = OrderedDict{String,Any}("name" => name, "kg" => kg, "mol" => mol)

"""
    feed_summary(feed; atomic_weights=nothing, species_elements=nothing) -> OrderedDict

Numbers for the live Feed summary of the Case tab:

    mode         "elements" or "species"
    rows         [{"name", "kg", "mol"}]: Total, Water, Solids, Ash, Ash-free solids
                 (species mode: Total, Water, Solids (not water), each species),
                 then each element; "mol" is `nothing` where it has no meaning
    elements     {element => {"kg", "mol"}} (without the water's H and O)
    used         {element => % used} (elemental mode)
    entered_sum  sum of the entered element percentages (elemental mode)

and, by name, total_kg, water_kg, water_mol, solids_kg, ash_kg, ash_free_kg.
Throws a FeedError naming the first problem.
"""
function feed_summary(feed, atomic_weights_=nothing, species_elements_=nothing;
                      atomic_weights=atomic_weights_, species_elements=species_elements_)
    aw = _atomic_weights(atomic_weights)
    M_water = water_molar_mass(aw)
    if feed_mode(feed) == "species"
        species_elements === nothing &&
            throw(FeedError("Species amounts need the species data (species_elements)."))
        molar_mass = molar_mass_function(species_elements, aw)
        masses = species_masses(feed[SPECIES_KEY], molar_mass)
        water_kg = masses["water"][1]
        dry_kg = sum(kg for (name, (kg, _)) in masses if name != "water"; init=0.0)
        element_kg = OrderedDict{String,Float64}()
        for (name, (kg, mol)) in masses
            name == "water" && continue
            for (element, count) in species_elements(name)
                element_kg[element] = get(element_kg, element, 0.0) + mol * count * aw[element] / 1000.0
            end
        end
        elements = OrderedDict{String,Any}(e => OrderedDict{String,Any}("kg" => kg, "mol" => 1000.0 * kg / aw[e])
                                           for (e, kg) in element_kg)
        rows = [_row("Total", water_kg + dry_kg, nothing),
                _row("Water", water_kg, 1000.0 * water_kg / M_water),
                _row("Solids (not water)", dry_kg, nothing)]
        append!(rows, [_row(name, kg, mol) for (name, (kg, mol)) in masses if name != "water"])
        append!(rows, [_row(e, v["kg"], v["mol"]) for (e, v) in elements])
        return OrderedDict{String,Any}("mode" => "species", "rows" => rows, "elements" => elements,
                                       "total_kg" => water_kg + dry_kg, "water_kg" => water_kg,
                                       "water_mol" => 1000.0 * water_kg / M_water, "solids_kg" => dry_kg,
                                       "ash_kg" => 0.0, "ash_free_kg" => dry_kg)
    end
    masses = elemental_masses(feed)
    elements_in = get(feed, "elements_wt_percent", nothing)
    used, entered_sum = scaled_elements(py_truth(elements_in) ? elements_in : _Dict())
    missing_ = [e for e in keys(used) if !haskey(aw, e)]
    isempty(missing_) ||
        throw(FeedError("Unknown element(s) $(join(missing_, ", ")); give their atomic weights."))
    elements = OrderedDict{String,Any}(
        e => OrderedDict{String,Any}("kg" => masses["ash_free"] * p / 100.0,
                                     "mol" => 1000.0 * masses["ash_free"] * p / 100.0 / aw[e])
        for (e, p) in used)
    rows = [_row("Total", masses["total"], nothing),
            _row("Water", masses["water"], 1000.0 * masses["water"] / M_water),
            _row("Solids", masses["solids"], nothing),
            _row("Ash", masses["ash"], nothing),
            _row("Ash-free solids", masses["ash_free"], nothing)]
    append!(rows, [_row(e, v["kg"], v["mol"]) for (e, v) in elements])
    return OrderedDict{String,Any}("mode" => "elements", "rows" => rows, "elements" => elements, "used" => used,
                                   "entered_sum" => entered_sum, "total_kg" => masses["total"],
                                   "water_kg" => masses["water"], "water_mol" => 1000.0 * masses["water"] / M_water,
                                   "solids_kg" => masses["solids"], "ash_kg" => masses["ash"],
                                   "ash_free_kg" => masses["ash_free"])
end

"""
    feed_fractions(feed_model, molar_mass=nothing) -> (s, a)

s = solids / total mass and a = ash / solids for a model Feed. Results use them
to convert per kg dry feed into per kg total feed (s) or per kg dry ash-free
feed (a). For an elemental analysis the ash is the dry mass outside the
elements (as in the legacy conversion). Species amounts are resolved with
`molar_mass(name) -> g/mol`; without it (or when a molar mass is unknown) s is
`nothing`.
"""
function feed_fractions(feed_model::Feed, molar_mass=nothing)
    if feed_model.species_amounts !== nothing
        molar_mass === nothing && return nothing, 0.0
        try
            feed_model = resolved(feed_model, molar_mass)
        catch e
            e isa Union{KeyError,ValueError,ArgumentError,MethodError,DivideError} || rethrow()
            return nothing, 0.0
        end
    end
    s = 1.0 / (1.0 + total_water_g_per_g_dry(feed_model))
    feed_model.species_moles === nothing || return s, 0.0
    if py_truth(feed_model.unmodelled_wt_percent)
        ash = feed_model.ash_wt_percent
        return s, (py_truth(ash) ? ash : 0.0) / 100.0
    end
    rest = 100.0 - sum(values(dry_element_wt_percent(feed_model, nothing)); init=0.0)
    return s, (rest > 0.0 ? rest : 0.0) / 100.0   # Python max(0.0, rest)
end

# ---------------------------------------------------------------- legacy feeds
const UNMODELLED_MESSAGE = "This case uses unmodelled_wt_percent, which the form can't show. " *
                           "Edit the file directly."
const SPECIES_NEEDS_DATA_MESSAGE = "This case gives the feed as species moles; converting it needs " *
                                   "the species data."

"""
1 + value, written without binary noise when value is a short decimal (Python:
float(Decimal(repr(float(value))) + 1), with Decimal's 28 significant digits).
"""
function _plus_one(value)
    x = py_float(value)
    isfinite(x) || return 1.0 + x
    m = match(r"^(-?)(\d+)(?:\.(\d*))?(?:e([+-]\d+))?$", py_repr(x))
    frac = something(m[3], "")
    coefficient = parse(BigInt, m[2] * frac) * (m[1] == "-" ? -1 : 1)
    exponent = (m[4] === nothing ? 0 : parse(Int, m[4])) - length(frac)
    if exponent >= 0
        c, e = coefficient * big(10)^exponent + 1, 0
    else
        c, e = coefficient + big(10)^(-exponent), exponent
    end
    digits = ndigits(abs(c))
    if digits > 28   # round to Decimal's default precision, half to even
        shift = digits - 28
        q, r = divrem(abs(c), big(10)^shift)
        half = big(10)^shift
        if 2r > half || (2r == half && isodd(q))
            q += 1
        end
        c, e = sign(c) * q, e + shift
    end
    return parse(Float64, string(c) * "e" * string(e))
end

"""
    legacy_feed_to_v2(feed, molar_mass=nothing) -> OrderedDict

Format-2 [feed] table for a legacy [feed] table (same model numbers). Throws a
`LegacyConversionError` when the feed cannot be written in format 2:
unmodelled_wt_percent > 0, or a species feed without `molar_mass`.

Elemental feeds use 1 kg of solids. The ash written is 100 minus the dry
element percentages, so that the dry elemental analysis the model uses is
unchanged: for a "dry" analysis without ash that sums to less than 100, the
rest becomes ash; where the old file gave an ash value that differs from that
rest (the old loader allowed 0.05 %), the rest is used.
"""
function legacy_feed_to_v2(feed, molar_mass=nothing)
    feed = deepcopy(_Dict(feed))
    salts = pop!(feed, "salts", nothing)
    py_truth(get(feed, "unmodelled_wt_percent", nothing)) && throw(LegacyConversionError(UNMODELLED_MESSAGE))
    model = _feed_from_kwargs(feed)
    water = total_water_g_per_g_dry(model)
    if model.species_moles !== nothing
        molar_mass === nothing && throw(LegacyConversionError(SPECIES_NEEDS_DATA_MESSAGE))
        dry_g = sum(py_float(n) * molar_mass(name) for (name, n) in model.species_moles; init=0.0)
        species = _Dict("water" => _Dict("value" => water * dry_g / 1000.0, "unit" => "kg"))
        for (name, n) in model.species_moles
            species[name] = _Dict("value" => py_float(n), "unit" => "mol")
        end
        out = _Dict("species" => species)
    else
        dry = dry_element_wt_percent(model, nothing)
        ash_given = model.ash_wt_percent
        divisor = 1.0
        if model.basis == "dry" && py_truth(ash_given)
            divisor = 1.0 - py_float(ash_given) / 100.0
        end
        elements = model.basis == "dry_ash_free" ? _Dict(model.elements_wt_percent) :
                   _Dict(e => (divisor != 1.0 ? py_float(v) / divisor : v) for (e, v) in model.elements_wt_percent)
        rest = 100.0 - sum(values(dry); init=0.0)
        out = _Dict("total_kg" => (!py_truth(model.moisture_wt_percent) ? _plus_one(water) : 1.0 + water),
                    "solids" => _Dict("value" => 1.0, "unit" => "kg"))
        if ash_given !== nothing && abs(rest - py_float(ash_given)) <= 1e-9 * max(1.0, abs(rest))
            out["ash"] = _Dict("value" => ash_given, "unit" => "%")
        elseif rest > 1e-9
            out["ash"] = _Dict("value" => parse(Float64, py_format_g(rest, 14)), "unit" => "%")
        end
        out["elements_wt_percent"] = elements
    end
    py_truth(salts) && (out["salts"] = salts)
    return out
end
