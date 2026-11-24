# Comparison of the case/feed/units port with the Python app: reads
# test/data/cases_python.json (written by test/dump_cases.py, not committed:
# ~0.5 MB) and builds the same dump on the Julia side (see dump_cases.py for
# what is in it). Included by test_units.jl, test_case.jl and test_feed.jl.
using HTGEquilibrium
using HTGEquilibrium: OrderedDict, ValueError, py_truth
using HTGEquilibrium.Units: UnitError
using JSON
using Test

const CASES_JSON = joinpath(@__DIR__, "data", "cases_python.json")
const PY_CASES = isfile(CASES_JSON) ? JSON.parsefile(CASES_JSON) : nothing
const REPO_ROOT = realpath(joinpath(@__DIR__, "..", "..", ".."))
PY_CASES === nothing && @info "test/data/cases_python.json not found: the comparisons with Python are skipped. " *
                              "Write it with `uv run python julia/HTGEquilibrium/test/dump_cases.py`."

# ---------------------------------------------------------------- the dump format
"""Julia value -> the JSON form dump_cases.py writes (dict order kept)."""
enc(::Nothing) = nothing
enc(x::Bool) = x
enc(x::Integer) = Int(x)
enc(x::AbstractFloat) = isfinite(x) ? Float64(x) : OrderedDict{String,Any}("__float__" => py_repr(x))
enc(x::AbstractString) = String(x)
enc(x::Union{AbstractVector,Tuple}) = Any[enc(v) for v in x]
function enc(x::AbstractDict)
    all(k -> k isa AbstractString, keys(x)) && return OrderedDict{String,Any}(k => enc(v) for (k, v) in x)
    return OrderedDict{String,Any}("__pairs__" => Any[Any[enc(k), enc(v)] for (k, v) in x])
end
enc(x::Union{Case,Feed,Salt}) = OrderedDict{String,Any}(string(f) => enc(getfield(x, f)) for f in fieldnames(typeof(x)))

"""JSON (as parsed) -> plain Julia data: OrderedDicts, Vector{Any}, Int, Float64."""
decode(x::AbstractDict) = haskey(x, "__float__") ? parse(Float64, replace(x["__float__"], "inf" => "Inf")) :
                          OrderedDict{String,Any}(String(k) => decode(v) for (k, v) in x)
decode(x::AbstractVector) = Any[decode(v) for v in x]
decode(x) = x

# Python exception type names of the Julia errors.
pytype(e::CaseError) = "ValueError"
pytype(e::UnitError) = "ValueError"
pytype(e::ValueError) = string(nameof(typeof(e)))   # FeedError, LegacyConversionError
pytype(e::ArgumentError) = "TypeError"   # Python TypeErrors (Feed(**data) with unknown keys, ...)
pytype(e) = "julia " * string(nameof(typeof(e)))

errdict(e) = OrderedDict{String,Any}("type" => pytype(e), "message" => e isa Union{ValueError,ArgumentError} ? e.msg : sprint(showerror, e))

"""{"ok": f(args...)} or {"error": {type, message}} (Python: attempt)."""
function attempt(f, args...; kwargs...)
    try
        return OrderedDict{String,Any}("ok" => f(args...; kwargs...))
    catch e
        e isa Union{ValueError,KeyError,ArgumentError,MethodError,DivideError} || rethrow()
        return OrderedDict{String,Any}("error" => errdict(e))
    end
end

# ---------------------------------------------------------------- comparison
# Python is the reference, not bit for bit (owner, PR #7): floats agree within
# FLOAT_RTOL relative, or FLOAT_ATOL absolute where Python's result is a
# cancellation residue (100 - sum of weight percents: 0.0 in Python, ~1e-14 in
# Julia); strings, ints, Int-vs-float, None and key order exactly.
const FLOAT_RTOL = 1e-12
const FLOAT_ATOL = 1e-13
mutable struct FloatStats
    compared::Int
    inexact::Int
    worst::Float64       # largest relative deviation among those within FLOAT_RTOL
    worst_abs::Float64   # largest absolute deviation among those only within FLOAT_ATOL
end
const FLOATS = FloatStats(0, 0, 0.0, 0.0)

"""Differences between a Julia dump `a` and the Python JSON `b` (dict key order
included). Floats within FLOAT_RTOL relative or FLOAT_ATOL absolute; the largest deviations are kept in FLOATS."""
function pydiff(a, b, path="", out=String[])
    if b isa AbstractDict
        if !(a isa AbstractDict)
            push!(out, "$path: Julia $(repr(a)), Python $(repr(b))")
        else
            ka, kb = collect(String.(keys(a))), collect(String.(keys(b)))
            ka == kb || push!(out, "$path: keys Julia $ka, Python $kb")
            for k in kb
                haskey(a, k) && pydiff(a[k], b[k], "$path.$k", out)
            end
        end
    elseif b isa AbstractVector
        if !(a isa AbstractVector) || length(a) != length(b)
            push!(out, "$path: Julia $(repr(a)), Python $(repr(b))")
        else
            for i in eachindex(b)
                pydiff(a[i], b[i], "$path[$i]", out)
            end
        end
    elseif b isa Bool
        a === b || push!(out, "$path: Julia $(repr(a)), Python $b")
    elseif b isa Integer
        (a isa Integer && !(a isa Bool) && a == b) || push!(out, "$path: Julia $(repr(a)), Python int $b")
    elseif b isa AbstractFloat
        if !(a isa AbstractFloat)
            push!(out, "$path: Julia $(repr(a)), Python float $b")
        else
            FLOATS.compared += 1
            if a !== b
                rel = abs(a - b) / max(abs(a), abs(b))
                FLOATS.inexact += 1
                if rel <= FLOAT_RTOL
                    FLOATS.worst = max(FLOATS.worst, rel)
                elseif abs(a - b) <= FLOAT_ATOL
                    FLOATS.worst_abs = max(FLOATS.worst_abs, abs(a - b))
                else
                    push!(out, "$path: Julia $a, Python $b")
                end
            end
        end
    elseif b isa AbstractString
        a == b || push!(out, "$path: Julia $(repr(a)), Python $(repr(b))")
    elseif b === nothing
        a === nothing || push!(out, "$path: Julia $(repr(a)), Python None")
    else
        push!(out, "$path: unexpected Python value $(repr(b))")
    end
    return out
end

macro test_same(julia, python, label="")
    quote
        local diffs = pydiff($(esc(julia)), $(esc(python)))
        isempty(diffs) || @error "Differences from Python $($(esc(label)))" diffs = first(diffs, 20)
        @test isempty(diffs)
    end
end

# ---------------------------------------------------------------- species data (registry stand-in)
const SPECIES_TABLE = PY_CASES === nothing ? OrderedDict{String,Any}() :
    OrderedDict{String,Any}(k => decode(v["elements"]) for (k, v) in PY_CASES["species"])
const SOLID_NAMES = PY_CASES === nothing ? Set{String}() :
    Set{String}(k for (k, v) in PY_CASES["species"] if v["solid"])
const SPECIES_ELEMENTS = species_elements_function(SPECIES_TABLE)

case_molar_mass(case=nothing) = molar_mass_function(SPECIES_ELEMENTS, case === nothing ? nothing : case.atomic_weights)

carbon_solids(case) = OrderedDict{String,Any}(
    name => py_truth(get(SPECIES_TABLE[name], "C", 0)) for name in get(case.phases, "solid", String[]) if name in SOLID_NAMES)

# ---------------------------------------------------------------- the dump of one case
function case_dump(case::Case)
    fields = enc(case)
    # Julia only ([aqueous] water_model); Python has no such field. Tested in test_case.jl.
    fields["water_model"] == "SUPCRT92" && delete!(fields, "water_model")
    delete!(fields, "legacy_notes")   # Julia only (dropped solver settings); tested in test_case.jl
    comparison = fields["comparison"]
    if comparison isa AbstractDict && isabspath(get(comparison, "file", ""))
        comparison["file"] = relpath(comparison["file"], REPO_ROOT)
    end
    out = OrderedDict{String,Any}("case" => fields)
    for (name, f) in (("temperatures", temperatures), ("pressures", pressures), ("grid", grid),
                      ("chains", chains), ("is_single_point", is_single_point))
        out[name] = enc(attempt(f, case))
    end
    ps = get(out["pressures"], "ok", nothing)
    ts = get(out["temperatures"], "ok", nothing)
    if py_truth(ts) && py_truth(ps)
        point = at(case, ts[1], ps[1])
        out["at"] = enc(OrderedDict(k => getfield(point, Symbol(k)) for k in ("temperature_C", "pressure_bar", "points", "refine")))
    end
    solids = carbon_solids(case)
    out["carbon_solids"] = enc(solids)
    out["carbon_bound_messages"] = Any[carbon_bound_messages(case), carbon_bound_messages(case, solids)]
    out["effective_carbon_bounds"] = enc([effective_carbon_bounds(case.carbon_bounds, true),
                                          effective_carbon_bounds(case.carbon_bounds, false)])
    out["fluid_phase_count"] = fluid_phase_count(case.phases, case.aqueous_enabled, case.aqueous_cutoff_C, ts)
    r = attempt(resolved, case.feed, case_molar_mass(case))
    numbers = OrderedDict{String,Any}("resolved" => enc(r))
    if haskey(r, "ok")
        feed = r["ok"]
        numbers["water"] = enc(attempt(total_water_g_per_g_dry, feed))
        feed.species_moles === nothing && (numbers["dry"] = enc(attempt(dry_element_wt_percent, feed, nothing)))
    end
    numbers["fractions"] = enc(attempt(feed_fractions, case.feed, case_molar_mass(case)))
    numbers["fractions_without_molar_mass"] = enc(attempt(feed_fractions, case.feed))
    out["feed"] = numbers
    if case.feed_v2 !== nothing
        out["feed_summary"] = enc(attempt(feed_summary, case.feed_v2, case.atomic_weights, SPECIES_ELEMENTS))
    end
    text = attempt(to_toml, case)
    haskey(text, "ok") && (text["ok"] = replace(text["ok"], REPO_ROOT => "<repo>"))
    out["to_toml"] = text
    return out
end

function converted_dump!(out, data, name)
    converted = attempt(legacy_to_v2, data, case_molar_mass())
    out["legacy_to_v2"] = enc(converted)
    # Python's dict == treats a NaN that is the same object as equal (deepcopy keeps floats): isequal
    if haskey(converted, "ok") && !isequal(converted["ok"], data)
        c = attempt(case_from_dict, converted["ok"]; name)
        out["converted"] = haskey(c, "ok") ? case_dump(c["ok"]) : c
    end
    return out
end

function text_dump(text, name)
    data = parse_toml(text)
    c = attempt(case_from_dict, data; name)
    out = haskey(c, "ok") ? case_dump(c["ok"]) : OrderedDict{String,Any}("error" => c["error"])
    return converted_dump!(out, data, name)
end

function file_dump(rel)
    path = joinpath(REPO_ROOT, rel)
    out = OrderedDict{String,Any}("path" => rel)
    c = attempt(load_case, path)
    haskey(c, "error") && (out["error"] = c["error"]; return out)
    merge!(out, case_dump(c["ok"]))
    return converted_dump!(out, read_toml(path), splitext(basename(path))[1])
end

function variant_text(base, replacements)
    text = read(joinpath(REPO_ROOT, base), String)
    for (old, new) in replacements
        @assert occursin(old, text) (base, old)
        text = replace(text, old => new)
    end
    return text
end

# The Python app still has the solvers "slsqp-fd" and "trust-constr" and the
# [solver] keys fd_ftol and fd_eps; the Julia version dropped them (owner, PR #10:
# an old case loads with "ipopt" and ignores the keys; see test_case.jl). The
# Python dump (dump_cases.py, unchanged) is mapped to what Julia gives instead:
# Case encodings lose the two fields (a dropped solver becomes "ipopt"), the
# to_toml texts lose their fd_ lines, and the solver-list message lists two
# solvers and the allowed [solver] keys without the two. Format-1 -> format-2 data
# (legacy_to_v2) is not mapped: it keeps the keys.
const PY_SOLVER_LIST = "solver must be one of ('ipopt', 'slsqp', 'trust-constr', 'slsqp-fd')"
const JL_SOLVER_LIST = "solver must be one of ('ipopt', 'slsqp')"
without_dropped_solvers(x) = x
without_dropped_solvers(x::AbstractVector) = Any[without_dropped_solvers(v) for v in x]
function without_dropped_solvers(x::AbstractString)
    x = replace(x, PY_SOLVER_LIST => JL_SOLVER_LIST)
    startswith(x, "Unknown [solver] key(s)") && (x = replace(x, "'fd_eps', 'fd_ftol', " => ""))
    occursin(r"(?m)^fd_(ftol|eps) = ", x) || return x
    x = replace(x, r"(?m)^fd_(ftol|eps) = [^\n]*\n" => "")
    return replace(x, r"\n\[solver\]\n(?=\n|$)" => "")   # a [solver] table left empty
end
function without_dropped_solvers(d::AbstractDict)
    out = OrderedDict{String,Any}(k => without_dropped_solvers(v) for (k, v) in d)
    if haskey(out, "temperature_C") && haskey(out, "fd_ftol")   # a Case encoding
        delete!(out, "fd_ftol")
        delete!(out, "fd_eps")
        get(out, "solver", nothing) in ("slsqp-fd", "trust-constr") && (out["solver"] = "ipopt")
    end
    return out
end

# Python's dump minus the keys only the Julia side does not write.
strip_keys(d::AbstractDict, ks...) = OrderedDict{String,Any}(k => v for (k, v) in d if !(k in ks))
