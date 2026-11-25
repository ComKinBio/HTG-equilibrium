# Lists of numbers typed in a text field (Python: htg/gui/editor.py parse_numbers,
# format_numbers, parse_points; spec WP1 item 4). In the core so that the app and
# the command line read lists the same way.
#
#     parse_numbers("300:400 n=5")            -> [300.0, 325.0, 350.0, 375.0, 400.0]
#     parse_numbers("24, 30"; unit="MPa")     -> [240.0, 300.0]   (converted to bar)
#     format_numbers([300, 350.5])            -> "300, 350.5"
#     parse_points("300, 240\n350 300")       -> [[300.0, 240.0], [350.0, 300.0]]

"""A list of numbers that cannot be read; the message can be shown under the field (Python: ValueError)."""
struct NumberListError <: ValueError
    msg::String
end

"""At most this many values from one range or count."""
const MAX_LIST_VALUES = 10000
const _NUMBER_FORMATS = "Write numbers like 300, 350, a range 300:400:50 or 300:400 n=5."

function _list_number(text)
    value = py_tryfloat(text)
    (value === nothing || !isfinite(value)) && throw(NumberListError("\"$text\" is not a number. $_NUMBER_FORMATS"))
    return value
end

"""Drop floating-point noise from generated or converted values: 299.99999999999994 -> 300."""
clean_number(value) = parse(Float64, py_format_g(value, 12))

# 'a:b:step' (count nothing) or 'a:b' with count -> values (in the unit typed)
function _list_range(token, count)
    parts = split(token, ':')
    if !(length(parts) in (2, 3)) || any(isempty, parts)
        throw(NumberListError("\"$token\" is not a range. Write start:stop:step (300:400:50) " *
                              "or start:stop n=count (300:400 n=5)."))
    end
    start, stop = _list_number(parts[1]), _list_number(parts[2])
    start > stop && throw(NumberListError("In \"$token\" the start is above the stop. Put the smaller value first."))
    if length(parts) == 3
        count === nothing || throw(NumberListError("\"$token n=$count\" gives both a step and a count. Use one of them."))
        step = _list_number(parts[3])
        step <= 0 && throw(NumberListError("In \"$token\" the step must be above zero."))
        steps = stop > start ? round(Int, (stop - start) / step) : 0
        steps + 1 > MAX_LIST_VALUES &&
            throw(NumberListError("\"$token\" gives more than $MAX_LIST_VALUES values. Use a larger step."))
        return [clean_number(start + i * step) for i in 0:steps if start + i * step <= stop + 1e-9 * step]
    end
    count === nothing && throw(NumberListError("\"$token\" needs a step or a count, e.g. $token:10 or $token n=5."))
    return vcat([clean_number(start + i * (stop - start) / (count - 1)) for i in 0:count-2], [stop])
end

function _list_count(token)
    text = token[3:end]
    if isempty(text) || !all(isdigit, text) || parse(BigInt, text) < 2
        throw(NumberListError("\"$token\": the count after n= must be a whole number of 2 or more."))
    end
    parse(BigInt, text) > MAX_LIST_VALUES &&
        throw(NumberListError("\"$token\": the count can be at most $MAX_LIST_VALUES."))
    return parse(Int, text)
end

"""
    parse_numbers(text; unit=nothing, sort=true) -> Vector{Float64}

Text of a number field -> numbers. Accepted (separators: comma, semicolon or
whitespace; formats can be mixed):

    300, 350, 400     list                               300, 350, 400
    300:400:50        range start:stop:step              300, 350, 400 (stops at or before stop)
    300:400 n=5       evenly spaced start:stop n=count   300, 325, 350, 375, 400
    250, 300:400:50   mixed                              250, 300, 350, 400

`n=` belongs to the `a:b` range just before it; the count is a whole number of
2 or more.

- `unit`: the unit the numbers are typed in (a unit of `Units`, e.g. "K" or
  "MPa"); the result is converted to the file unit (°C or bar). `nothing`: the
  numbers are already in file units.
- `sort`: sort and remove duplicates (temperature and pressure grids). false
  keeps the order as typed (e.g. one "T, P" point).

Throws a `NumberListError` (a `ValueError`) with a message that can be shown
under the field; an unknown unit is a `Units.UnitError`.
"""
function parse_numbers(text; unit=nothing, sort::Bool=true)
    text = replace(strip(something(text, "")), r"\s*:\s*" => ":")
    text = replace(text, r"\bn\s*=\s*" => "n=")
    isempty(text) && return Float64[]
    tokens = [String(t) for t in split(text, r"[,;\s]+") if !isempty(t)]
    out = Float64[]
    for (i, token) in enumerate(tokens)
        if startswith(token, "n=")
            previous = i > 1 ? tokens[i-1] : ""
            count(==(':'), previous) == 1 ||
                throw(NumberListError("\"$token\" must follow a range start:stop, e.g. 300:400 n=5."))
            continue
        end
        if occursin(':', token)
            following = i < length(tokens) ? tokens[i+1] : ""
            append!(out, _list_range(token, startswith(following, "n=") ? _list_count(following) : nothing))
        else
            push!(out, _list_number(token))
        end
    end
    if unit !== nothing && !(Units.unit_key(unit) in values(Units.FILE_UNITS))
        out = [clean_number(Units.to_file(v, unit)) for v in out]
    end
    return sort ? _sorted_unique(out) : out
end

# sorted(set(values)) of Python (0.0 and -0.0 are one value; the first is kept)
function _sorted_unique(values)
    out, seen = Float64[], Set{Float64}()
    for v in values
        key = v == 0 ? 0.0 : v
        key in seen && continue
        push!(seen, key)
        push!(out, v)
    end
    return Base.sort!(out)
end

"""
    format_numbers(value, unit=nothing) -> String

Stored value(s) in °C or bar -> field text in `unit` (`nothing`: as stored).
`nothing` gives "".
"""
function format_numbers(value, unit=nothing)
    value === nothing && return ""
    list = value isa AbstractVector ? value : [value]
    unit === nothing && return Units.format_list(list)
    return join((Units.display(v, unit) for v in list), ", ")
end

"""
    parse_points(text, temperature_unit=nothing, pressure_unit=nothing) -> Vector{Vector{Float64}}

Lines "T, P" (in the given units) -> [[T_C, P_bar], ...]; empty lines are skipped.
"""
function parse_points(text, temperature_unit=nothing, pressure_unit=nothing)
    points = Vector{Float64}[]
    for line in split(something(text, ""), r"\r\n|\r|\n")
        isempty(strip(line)) && continue
        numbers = parse_numbers(line; sort=false)
        length(numbers) == 2 || throw(NumberListError("\"$(strip(line))\": give a temperature and a pressure."))
        T, P = numbers
        temperature_unit === nothing || (T = clean_number(Units.to_file(T, temperature_unit)))
        pressure_unit === nothing || (P = clean_number(Units.to_file(P, pressure_unit)))
        push!(points, [T, P])
    end
    return points
end
