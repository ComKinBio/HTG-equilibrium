# The Case tab without the page: wording, field parsing, notices and derived
# numbers (Python: gui/case_form.py, and default_feed and scalar of gui/editor.py).
#
# Everything the Case tab shows that is computed (point count, the feed numbers,
# the carbon conflict notice, the Initial guesses tooltip, the SUPCRT92 band
# notice, the first problem for the status line) comes from the functions here,
# so it can be tested without a browser. The number lists (parse_numbers,
# format_numbers, parse_points) are in the core.

const Units = HTGEquilibrium.Units

# ---------------------------------------------------------------- wording (spec section 5.1)
const FEED_CAPTION = "The wet feed going into the reactor: solids plus water."
const FEED_CAPTION_SPECIES = "The wet feed going into the reactor: species from the database plus water. " *
                             "No ash in this mode."
const PRODUCTS_CAPTION = "Species that may form in each phase. Type to search. Ions can only go in the " *
                         "aqueous phase."
const CARBON_CAPTION = "Share of the feed carbon that ends up in each phase. Free = the solver decides."
const CUTOFF_CAPTION = "This is a modelling choice. If you know that the aqueous phase is definitely not " *
                       "present above a certain temperature, then marking it here helps the equilibrium " *
                       "solver. It is not necessary to input a value."
# Owner decision (PR #10): trust-constr and "SLSQP, finite differences" are gone, and so is their sentence.
const METHOD_TOOLTIP = "IPOPT is the default and usually the fastest. Try SLSQP if IPOPT fails at some points."
const CONTINUATION_TOOLTIP = "Start each point from the solution at the previous temperature (same " *
                             "pressure), or the previous point in a list. Usually faster and smoother " *
                             "along a sweep. Turn off to start every point fresh; then all points can " *
                             "run in parallel."
const ELEMENTS_CAPTION = "Scaled to sum to 100 % if needed"
const SPECIES_FOOTER = "Only species with complete data are listed"
const ADVANCED_CAPTION = "Leave a field empty to use the default"
const INITIAL_COMPOSITION_CAPTION = "Automatic picks liquid below the aqueous cut-off and vapour above it."
const SEED_CAPTION = "Set a number to get the same starting points every run."
const ENERGY_SCALE_CAPTION = "Scales the Gibbs energy for the optimizer. Change only if the solver struggles."
const DISTINCT_CAPTION = "Two results from different starting points count as different if either limit below is exceeded."
const COMPOSITION_TOL_CAPTION = "Different if any species amount differs by more than this."
const GIBBS_TOL_CAPTION = "Different if the Gibbs energies differ by more than this."
const ATOMIC_WEIGHTS_CAPTION = "Override the standard values [g / mol]. Affects the feed and element balances only."
const MEASURED_CAPTION = "The file needs a temperature column (T_C), then one column per species named " *
                         "phase:species, e.g. gas:methane."
const CHECK_CAPTION = "The checks need no equation of state. Passing does not guarantee that the solver converges."

# Water model (owner decision PR #7; Julia only). Wording proposed, see the README.
const WATER_MODEL_LABEL = "Water model"
const WATER_MODEL_TOOLTIP = "The water properties of the Electrolyte model: the ion data (HKF) and the " *
                            "dielectric constant. SUPCRT92 is the reference and the default. IAPWS-95 is the " *
                            "international standard for water and has no failure bands near the critical point."
const WATER_MODEL_MOLECULAR = "Used only with the Electrolyte model."
const WATER_MODEL_NO_AQUEOUS = "Used only when the case includes an aqueous phase."

const INPUT_FORMATS = [  # Format, Type, Gives
    ("List", "300, 350, 400", "300, 350, 400"),
    ("Range (start:stop:step)", "300:400:50", "300, 350, 400"),
    ("Evenly spaced (start:stop n=count)", "300:400 n=5", "300, 325, 350, 375, 400"),
    ("Mixed", "250, 300:400:50", "250, 300, 350, 400"),
]

const PHASE_NAMES = HTGEquilibrium.PHASE_LABELS
const CARBON_PHASES = ("gas", "solid", "aqueous", "organic_liquid")
const CARBON_CHOICES = OrderedDict("free" => "free", "exact" => "exactly", "min" => "at least", "max" => "at most")
const CARBON_RANGE = "range"   # min and max together: only offered when the file has it
const CARBON_RANGE_LABEL = "between"
const SOLID_CARBON = OrderedDict("char" => "char (solid carbon) only", "all_solids" => "all solids")
const AQUEOUS_MODELS = OrderedDict("electrolyte" => "Electrolyte (ePC-SAFT, ions allowed)",
                                   "molecular" => "Molecular (PC-SAFT, no ions)")
const WATER_MODEL_CHOICES = OrderedDict("SUPCRT92" => "SUPCRT92", "IAPWS-95" => "IAPWS-95")
# Owner decision (PR #10): only these two.
const METHODS = OrderedDict("ipopt" => "IPOPT (default)", "slsqp" => "SLSQP")
const INITIAL_COMPOSITIONS = OrderedDict("auto" => "Automatic", "liquid" => "Water mostly liquid (aqueous)",
                                         "vapor" => "Water mostly vapour (gas)", "lp" => "From a linear program")
const STABILITY_TEST = OrderedDict("default" => "On when using several initial guesses", "on" => "Always on",
                                   "off" => "Off")
const QUANTITY_CHOICES = OrderedDict("mol_per_kg_dry_feed" => "mol / kg dry feed",
                                     "mol_per_mol_C" => "mol / mol dry-feed C", "g_per_kg_dry_feed" => "g / kg dry feed")
const MODE_TOGGLE = OrderedDict("grid" => "Every T at every P", "points" => "List of (T, P) points")
const FEED_TOGGLE = OrderedDict("elements" => "Elemental analysis", "species" => "Species amounts")

const PERIODIC_TABLE = split(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se " *
    "Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb " *
    "Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm " *
    "Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og")

"""The species data are not loaded yet (the window opens first; see the README)."""
const LOADING_MESSAGE = "The app is still loading the species data. Try again in a moment."

# ---------------------------------------------------------------- case data
"""`data[path...]`, or `default` when a key is missing (Python: case_form.get)."""
function getpath(data, path, default=nothing)
    for key in path
        (data isa AbstractDict && haskey(data, key)) || return default
        data = data[key]
    end
    return data
end

"""Set `data[path...] = value`; `nothing` (or "") removes the key (Python: case_form.put)."""
function putpath!(data, path, value)
    for key in path[1:end-1]
        data = get!(data, key) do
            OrderedDict{String,Any}()
        end
    end
    if value === nothing || value == ""
        delete!(data, path[end])
    else
        data[path[end]] = value
    end
    return nothing
end

"""Drop floating-point noise from converted values: 299.99999999999994 -> 300."""
tidy(value, digits=12) = parse(Float64, HTGEquilibrium.py_format_g(value, digits))

# ---------------------------------------------------------------- one number in a text field
"""A field's text that cannot be used; the message goes under the field (Python: ValueError)."""
struct FieldError <: HTGEquilibrium.ValueError
    msg::String
end

"""
    parse_number(text, label=nothing; required=false, minimum=nothing, maximum=nothing, above=nothing,
                 integer=false) -> Float64, Int or nothing

Text of a number field -> a number (an Int if `integer`), or `nothing` when
empty and not required. `minimum`/`maximum` are inclusive limits, `above` an
exclusive lower limit. Throws a `FieldError` with a message for under the
field; with `label` the message starts with the field name (for the status line).
"""
function parse_number(text, label=nothing; required=false, minimum=nothing, maximum=nothing, above=nothing,
                      integer=false)
    fail(message) = throw(FieldError(label === nothing ? message : problem_text(label, message)))
    text = strip(text === nothing ? "" : string(text))
    if isempty(text)
        required && fail("Give a value.")
        return nothing
    end
    value = something(HTGEquilibrium.py_tryfloat(text), NaN)
    isfinite(value) || fail("\"$text\" is not a number.")
    if integer
        isinteger(value) || fail("Give a whole number.")
        value = Int(value)
    end
    minimum !== nothing && value < minimum && fail("Give a value of $(Units.format_number(minimum)) or more.")
    maximum !== nothing && value > maximum && fail("Give a value of $(Units.format_number(maximum)) or less.")
    above !== nothing && value <= above && fail("Give a value above $(Units.format_number(above)).")
    return value
end

"""A stored number as field text ("" for `nothing`)."""
number_text(value, unit=nothing) = value === nothing ? "" : Units.format_number(value, unit)

"""\"3 · Ash\", \"Ash can't be ...\" -> the message alone; otherwise \"Label: message\"."""
function problem_text(label, message)
    name = strip(last(split(first(split(label, " [")), " · ")))
    (isempty(name) || startswith(lowercase(message), lowercase(name))) && return message
    if length(message) >= 2 && isuppercase(message[1]) && islowercase(message[nextind(message, 1)])
        message = lowercase(message[1:1]) * message[nextind(message, 1):end]
    end
    return "$name: $message"
end

# ---------------------------------------------------------------- temperature and pressure
_values(value) = value === nothing ? Any[] : value isa AbstractVector ? collect(value) : Any[value]

"""Temperatures [°C] typed in [conditions] (grid or points), without the cut-off."""
function condition_temperatures(conditions)
    conditions = something(conditions, Dict())
    points = get(conditions, "points", nothing)
    points === nothing || return [Float64(p[1]) for p in points]
    return [Float64(t) for t in _values(get(conditions, "temperature_C", nothing))]
end

"""A Case of the conditions alone (empty feed), for counting points; `nothing` when incomplete."""
function conditions_case(conditions, cutoff_C=nothing)
    conditions = something(conditions, Dict())
    try
        return Case(; temperature_C=get(conditions, "temperature_C", nothing),
                    pressure_bar=get(conditions, "pressure_bar", nothing), points=get(conditions, "points", nothing),
                    feed=Feed(), phases=OrderedDict{String,Any}(), aqueous_cutoff_C=cutoff_C)
    catch exc
        exc isa Union{HTGEquilibrium.ValueError,ArgumentError,MethodError,KeyError,BoundsError,InexactError} || rethrow()
        return nothing
    end
end

"""
    point_count(conditions, cutoff_C=nothing) -> Int or nothing

Number of (T, P) points the case will solve, or `nothing` when the conditions
are incomplete. Counted as the calculation does (a cut-off inside the
temperature range adds a point).
"""
function point_count(conditions, cutoff_C=nothing)
    case = conditions_case(conditions, cutoff_C)
    case === nothing && return nothing
    try
        return length(grid(case))
    catch exc
        exc isa Union{HTGEquilibrium.ValueError,ArgumentError,MethodError,KeyError,BoundsError,InexactError} || rethrow()
        return nothing
    end
end

plural(count, word, words=nothing) = "$count $(count == 1 ? word : something(words, word * "s"))"

"""\"13 points\", \"1 point\"; \"\" when unknown."""
point_count_text(count) = count === nothing ? "" : plural(count, "point")

ready_text(points, pieces) = "Ready to run: $(plural(points, "point")) in $(plural(pieces, "independent piece"))."

"""The live tooltip of the Initial guesses field."""
function initial_guesses_tooltip(count)
    guesses = count == 1 ? "initial guess" : "initial guesses"
    return "You will start from $count $guesses. The optimization algorithm finds local equilibrium " *
           "points, but these may not be the global equilibrium. More starting points increase " *
           "confidence in finding the global equilibrium, at a higher computational cost. If several " *
           "distinct local equilibrium points are found, they are listed on the Results tab so you " *
           "can compare them."
end

_number_value(x) = x isa Real && !(x isa Bool)

"""
    conflict_notice(data, temperature_unit="C") -> String or nothing

The Aqueous phase card's conflict notice. Shown when (1) the cut-off has a
value, (2) the aqueous carbon choice is exactly or at least with a value above
0, and (3) at least one temperature of the case is at or above the cut-off.
"""
function conflict_notice(data, temperature_unit="C")
    cutoff = getpath(data, ("aqueous", "cutoff_C"))
    bound = something(getpath(data, ("constraints", "carbon", "aqueous")), Dict())
    (cutoff === nothing || !(bound isa AbstractDict)) && return nothing
    kind = haskey(bound, "exact") ? "exact" : haskey(bound, "min") ? "min" : nothing
    (kind === nothing || !_number_value(bound[kind]) || bound[kind] <= 0) && return nothing
    temperatures = try
        _number_value(cutoff) || return nothing
        condition_temperatures(get(data, "conditions", nothing))
    catch exc
        exc isa Union{MethodError,ArgumentError,InexactError,BoundsError} || rethrow()
        return nothing
    end
    any(t -> t >= cutoff, temperatures) || return nothing
    choice = kind == "exact" ? "exactly" : "at least"
    percent = Units.format_number(Units.fraction_to_percent(bound[kind]))
    T = "$(Units.display(Float64(cutoff), temperature_unit)) $(Units.label(temperature_unit))"
    return "Aqueous carbon is set to $choice $percent % (Carbon split). At $T and above there is no " *
           "aqueous phase, so that setting is ignored there. Below $T it applies."
end

# ---------------------------------------------------------------- SUPCRT92 failure bands (owner decision PR #4)
"""
    band_points(data) -> Vector{Tuple{Float64,Float64}}

The (T °C, P bar) points of the case where the electrolyte phase uses SUPCRT92:
none unless the aqueous phase is included, the model is Electrolyte and the
water model SUPCRT92; points at or above the cut-off have no aqueous phase and
are left out. Empty when the conditions are incomplete.
"""
function band_points(data)
    getpath(data, ("aqueous", "enabled"), true) === false && return Tuple{Float64,Float64}[]
    getpath(data, ("aqueous", "model"), "electrolyte") == "electrolyte" || return Tuple{Float64,Float64}[]
    getpath(data, ("aqueous", "water_model"), "SUPCRT92") == "SUPCRT92" || return Tuple{Float64,Float64}[]
    cutoff = getpath(data, ("aqueous", "cutoff_C"))
    case = conditions_case(get(data, "conditions", nothing))
    case === nothing && return Tuple{Float64,Float64}[]
    points = try
        grid(case)
    catch exc
        exc isa Union{HTGEquilibrium.ValueError,ArgumentError,MethodError,KeyError,BoundsError,InexactError} || rethrow()
        return Tuple{Float64,Float64}[]
    end
    return [(Float64(T), Float64(P)) for (T, P) in points if !(_number_value(cutoff) && T >= cutoff)]
end

"""
    band_hits(points; check=supcrt92_band_check, margin_C=0.5) -> Vector{Tuple{Float64,Float64}}

The points that are in or within `margin_C` of a SUPCRT92 failure band at
their pressure (`supcrt92_band_check`, about 0.3 s for each new pressure).
"""
function band_hits(points; check=supcrt92_band_check, margin_C=0.5)
    isempty(points) && return Tuple{Float64,Float64}[]
    result = check([T for (T, _) in points], [P for (_, P) in points]; margin_C)
    return [(r.T_C, r.P_bar) for r in result if r.status != :ok]
end

const BAND_LIST_MAX = 8

"""
    band_notice(hits, temperature_unit="C", pressure_unit="bar") -> String or nothing

The notice of the Temperature and pressure card for the points in or near a
SUPCRT92 failure band (warn only: the temperatures are not moved or dropped).
Wording proposed, see the README.
"""
function band_notice(hits, temperature_unit="C", pressure_unit="bar")
    isempty(hits) && return nothing
    groups = OrderedDict{Float64,Vector{Float64}}()
    for (T, P) in sort(hits; by=h -> (h[2], h[1]))
        push!(get!(groups, P, Float64[]), T)
    end
    parts = String[]
    shown = 0
    for (P, temps) in groups
        shown >= BAND_LIST_MAX && break
        temps = temps[1:min(end, BAND_LIST_MAX - shown)]
        shown += length(temps)
        push!(parts, join((Units.display(T, temperature_unit) for T in temps), ", ") *
                     " $(Units.label(temperature_unit)) at $(Units.display(P, pressure_unit)) $(Units.label(pressure_unit))")
    end
    more = length(hits) - shown
    list = join(parts, "; ") * (more > 0 ? " and $more more" : "")
    return "SUPCRT92 can fail at $list: in or within 0.5 $(Units.label(temperature_unit)) of a band where its " *
           "water density does not converge. Results at these points may be wrong. Move these temperatures, " *
           "or choose IAPWS-95 as the water model."
end

# ---------------------------------------------------------------- carbon split
"""
    carbon_choice(bound) -> (choice, first fraction, second fraction)

`choice` is a key of CARBON_CHOICES, or CARBON_RANGE for min and max together
(first = min, second = max).
"""
function carbon_choice(bound)
    (!(bound isa AbstractDict) || isempty(bound)) && return "free", nothing, nothing
    haskey(bound, "exact") && return "exact", bound["exact"], nothing
    haskey(bound, "min") && haskey(bound, "max") && return CARBON_RANGE, bound["min"], bound["max"]
    haskey(bound, "min") && return "min", bound["min"], nothing
    return "max", get(bound, "max", nothing), nothing
end

"""The [constraints.carbon] entry for a row of the Carbon split card (`nothing` = free)."""
function carbon_bound(choice, first=nothing, second=nothing)
    choice == "free" && return nothing
    choice == CARBON_RANGE && return OrderedDict{String,Any}("min" => first, "max" => second)
    return OrderedDict{String,Any}(choice => first)
end

# ---------------------------------------------------------------- feed
"""\"c\" -> \"C\", \" na \" -> \"Na\"."""
function element_symbol(text)
    text = strip(text === nothing ? "" : string(text))
    isempty(text) && return ""
    return uppercase(text[1:1]) * lowercase(text[nextind(text, 1):end])
end

"""Why `symbol` cannot be added to an element table (`nothing` when it can)."""
function element_error(symbol, existing=())
    isempty(symbol) && return "Type an element symbol, e.g. Cl."
    symbol in PERIODIC_TABLE || return "\"$symbol\" is not an element symbol."
    symbol in existing && return "$symbol is already in the table."
    return nothing
end

"""The scaling notice, or `nothing` when the entered percentages add up to 100.00."""
function entered_sum_notice(entered_sum)
    text = @sprintf("%.2f", entered_sum)
    text == "100.00" && return nothing
    return "The element percentages add up to $text %. They have been scaled to 100 % (the Used column)."
end

"""The field a FeedError message is about: "total", "solids", "ash", "elements",
"element:X", "species", "species:X" (or "feed")."""
function feed_error_field(message)
    for (prefix, field) in (("Total mass", "total"), ("Solids", "solids"), ("Ash", "ash"),
                            ("Elements", "elements"), ("Unknown element", "elements"))
        startswith(message, prefix) && return field
    end
    startswith(message, "Element ") && return "element:" * first(split(message[length("Element ")+1:end], ':'))
    startswith(message, "Species amounts: water") && return "species:water"
    startswith(message, "Species amounts") && return "species"
    startswith(message, "Species ") && return "species:" * first(split(message[length("Species ")+1:end], ':'))
    if startswith(message, "Unknown feed species")
        return occursin('\'', message) ? "species:" * split(message, '\'')[2] : "species"
    end
    if occursin("(species ", message)
        return "species:" * rstrip(last(split(message, "(species ")), ['.', ')'])
    end
    return "feed"
end

"""
    check_feed(feed, atomic_weights=nothing, species_elements=nothing) -> (summary, field, message)

The feed summary, or `nothing` with the first problem.
"""
function check_feed(feed, atomic_weights=nothing, species_elements=nothing)
    try
        return feed_summary(something(feed, OrderedDict{String,Any}()), atomic_weights, species_elements), nothing, nothing
    catch exc
        exc isa FeedError || rethrow()
        return nothing, feed_error_field(exc.msg), exc.msg
    end
end

"""
A "%" or "kg" amount in the other unit; `whole_kg` is what the % refers to.
Returns the value unchanged when it cannot be converted (no value or no whole).
"""
function convert_share(value, old_unit, new_unit, whole_kg)
    (value === nothing || old_unit == new_unit || whole_kg === nothing || whole_kg == 0) && return value
    converted = new_unit == "kg" ? value / 100.0 * whole_kg : 100.0 * value / whole_kg
    return tidy(converted, 6)
end

"""A species amount in "kg" or "mol" in the other unit (`molar_mass` in g/mol)."""
function convert_species_amount(value, old_unit, new_unit, molar_mass)
    (value === nothing || old_unit == new_unit || molar_mass === nothing || molar_mass == 0) && return value
    converted = new_unit == "mol" ? value * 1000.0 / molar_mass : value * molar_mass / 1000.0
    return tidy(converted, 6)
end

"""Feed summary masses: 8.500, 0.127; small values with 3 significant figures."""
function format_kg(value)
    value == 0 && return "0"
    abs(value) >= 0.01 && return @sprintf("%.3f", value)
    return Units.format_number(value, nothing, 3)
end

"""Feed summary amounts: 471.8, 37.96 (4 significant figures, trailing zeros kept)."""
function format_mol(value)
    value === nothing && return ""
    value == 0 && return "0"
    abs(value) >= 1000 && return @sprintf("%.0f", value)
    abs(value) < 1e-3 && return Units.format_number(value, nothing, 4)
    # Python f"{value:#.4g}".rstrip("."): here the exponent is between -3 and 3
    e = @sprintf("%.3e", value)
    exponent = parse(Int, e[findlast('e', e)+1:end])
    exponent >= 4 && return @sprintf("%.0f", value)   # 999.95 rounds up to 1000
    return String(rstrip(Printf.format(Printf.Format("%.$(3 - exponent)f"), value), '.'))
end

"""Used column and sums: two decimals."""
format_percent(value) = @sprintf("%.2f", value)

# ---------------------------------------------------------------- file bar
"""\"cases/manure.toml · unsaved changes\"."""
function file_status(path, dirty, cwd=nothing)
    path === nothing && return "New case · not saved yet"
    path = resolve_path(path)
    base = resolve_path(something(cwd, pwd()))
    rel = relpath(path, base)
    where = startswith(rel, "..") || isabspath(rel) ? path : rel
    return dirty ? "$where · unsaved changes" : where
end

# ---------------------------------------------------------------- from gui/editor.py
"""A list of one value -> the value (one temperature or pressure is written as a number)."""
scalar(value) = value isa AbstractVector && length(value) == 1 ? value[1] : value

"""A starting feed for a mode (the template's values)."""
function default_feed(mode, feed_species=())
    if mode == "species"
        species = OrderedDict{String,Any}("water" => OrderedDict{String,Any}("value" => 4.5, "unit" => "kg"))
        "cellulose" in feed_species && (species["cellulose"] = OrderedDict{String,Any}("value" => 0.5, "unit" => "kg"))
        return OrderedDict{String,Any}("species" => species)
    end
    return OrderedDict{String,Any}(
        "total_kg" => 5.0, "solids" => OrderedDict{String,Any}("value" => 20.0, "unit" => "%"),
        "ash" => OrderedDict{String,Any}("value" => 8.44, "unit" => "%"),
        "elements_wt_percent" => OrderedDict{String,Any}("C" => 30.4, "H" => 4.41, "O" => 53.58, "N" => 2.66,
                                                         "S" => 0.51))
end
