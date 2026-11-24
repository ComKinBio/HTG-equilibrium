# Small pieces of Python behaviour the port depends on, so that Julia reads the
# same files into the same values, with the same error messages:
#
#   ValueError            abstract type: Python ValueErrors of the case/feed/units port
#   py_repr, py_str       repr() / str() (floats: shortest round trip in Python's layout)
#   py_float, py_tryfloat Python float(x) (numbers, bools, strings in Python's grammar)
#   py_int                Python int(x)
#   py_round              round(x, n) of a float (correctly rounded decimal)
#   py_format_g           f"{x:.{p}g}" (C's %g; Julia's Printf %g differs)
#   py_isclose            math.isclose
#   py_truth              truthiness (None, 0, "", empty containers are false)
#   py_sorted             sorted() of strings
#   pypow                 x ** y for floats (C pow)
#   isna                  pandas isna
#   casefold              str.casefold()
#   py_slice              s[a:b] on characters (0-based, clamped)
#   read_text             open(path, encoding='utf-8', errors='ignore').read() (universal newlines)
#   py_splitlines         str.splitlines(keepends)
#
# TOML in file order (Python's tomllib) is in toml_ordered.jl (`parse_toml`, `read_toml`).
#
# Error conventions: the case/feed/units port throws subtypes of `ValueError`
# (CaseError, FeedError, LegacyConversionError, Units.UnitError) where Python
# raises ValueError. The data layer (databases, registry, user species) throws
# ArgumentError with Python's ValueError/TypeError text, KeyError(message) for
# KeyError. Python TypeErrors are ArgumentErrors everywhere.

"""
Python's `ValueError` family. Error types that are ValueErrors in Python
(`CaseError`, `FeedError`, `LegacyConversionError`, `Units.UnitError`)
subtype it, so `e isa ValueError` catches what `except ValueError` catches
in the case/feed code. (The data layer still uses ArgumentError; see README.)
"""
abstract type ValueError <: Exception end
Base.showerror(io::IO, e::ValueError) = print(io, e.msg)

"""Python `x ** y` for floats: C `pow` (the platform libm, as CPython)."""
pypow(x::Real, y::Real) = ccall(:pow, Float64, (Float64, Float64), Float64(x), Float64(y))

isna(x) = x === missing || x === nothing || (x isa AbstractFloat && isnan(x))   # pd.isna

casefold(s::AbstractString) = Base.Unicode.normalize(String(s); casefold=true, compose=false)

# ---------------------------------------------------------------- repr, str
"""
    py_repr(x) -> String

Python's `repr` for the values case files and databases hold: None, bools,
ints, floats (shortest round-trip digits, Python's switch to exponent notation:
`1e-05`, `1e+16`, `100.0`), strings (Python's quoting), lists/vectors, tuples
and dicts.
"""
py_repr(::Nothing) = "None"
py_repr(x::Bool) = x ? "True" : "False"
py_repr(x::Integer) = string(x)
py_repr(x::AbstractFloat) = _pyfloat_repr(Float64(x))
py_repr(x::AbstractString) = _pystr_repr(String(x))
py_repr(x::Symbol) = _pystr_repr(String(x))
py_repr(x::AbstractVector) = "[" * join((py_repr(v) for v in x), ", ") * "]"
py_repr(x::Tuple) = length(x) == 1 ? "(" * py_repr(x[1]) * ",)" : "(" * join((py_repr(v) for v in x), ", ") * ")"
py_repr(x::AbstractDict) = "{" * join((py_repr(k) * ": " * py_repr(v) for (k, v) in x), ", ") * "}"
py_repr(x) = repr(x)

"""Python's `str()`: like `py_repr`, except that strings are not quoted (and pandas NaN is "nan")."""
py_str(x::AbstractString) = String(x)
py_str(::Missing) = "nan"
py_str(x::Union{Nothing,Bool,Integer,AbstractFloat,AbstractVector,Tuple,AbstractDict}) = py_repr(x)
py_str(x) = string(x)

function _pyfloat_repr(x::Float64)
    isnan(x) && return "nan"
    isinf(x) && return x > 0 ? "inf" : "-inf"
    x == 0 && return signbit(x) ? "-0.0" : "0.0"
    m = match(r"^(\d+)(?:\.(\d+))?(?:e([+-]?\d+))?$", string(abs(x)))  # shortest digits
    frac = something(m[2], "")
    digits = m[1] * frac
    exp10 = (m[3] === nothing ? 0 : parse(Int, m[3])) - length(frac)  # x = digits * 10^exp10
    digits = lstrip(digits, '0')
    stripped = rstrip(digits, '0')
    exp10 += length(digits) - length(stripped)
    digits = String(stripped)
    decpt = length(digits) + exp10          # x = 0.digits * 10^decpt
    n = length(digits)
    text = if -4 < decpt <= 16
        if decpt <= 0
            "0." * "0"^(-decpt) * digits
        elseif decpt >= n
            digits * "0"^(decpt - n) * ".0"
        else
            digits[1:decpt] * "." * digits[decpt+1:end]
        end
    else
        e = decpt - 1
        (n > 1 ? digits[1:1] * "." * digits[2:end] : digits) * "e" * (e < 0 ? "-" : "+") * lpad(abs(e), 2, '0')
    end
    return signbit(x) ? "-" * text : text
end

function _pystr_repr(s::String)
    quote_char = ('\'' in s && !('"' in s)) ? '"' : '\''
    io = IOBuffer()
    print(io, quote_char)
    for c in s
        if c == quote_char || c == '\\'
            print(io, '\\', c)
        elseif c == '\n'
            print(io, "\\n")
        elseif c == '\r'
            print(io, "\\r")
        elseif c == '\t'
            print(io, "\\t")
        elseif !isprint(c)
            u = UInt32(c)
            print(io, u <= 0xff ? "\\x" * string(u, base=16, pad=2) :
                      u <= 0xffff ? "\\u" * string(u, base=16, pad=4) : "\\U" * string(u, base=16, pad=8))
        else
            print(io, c)
        end
    end
    print(io, quote_char)
    return String(take!(io))
end

# ---------------------------------------------------------------- format
"""
    py_format_g(value, p) -> String

Python's `f"{value:.{p}g}"` (C's "%.<p>g"). Julia's Printf `%g` picks fixed or
exponent notation before rounding (999999.5 with p = 6 gives "1000000"; C and
Python give "1e+06"), so the exponent is taken from the rounded "%e" form, as the
C standard specifies. Trailing zeros are removed.
"""
function py_format_g(value::Real, p::Integer)
    value = Float64(value)
    isfinite(value) || return isnan(value) ? "nan" : value > 0 ? "inf" : "-inf"
    p = max(p, 1)
    e = Printf.format(Printf.Format("%.$(p - 1)e"), value)
    x = parse(Int, e[findlast('e', e)+1:end])
    if -4 <= x < p
        text = Printf.format(Printf.Format("%.$(p - 1 - x)f"), value)
        return occursin('.', text) ? String(rstrip(rstrip(text, '0'), '.')) : text
    end
    mantissa, exponent = split(e, 'e')
    mantissa = occursin('.', mantissa) ? rstrip(rstrip(mantissa, '0'), '.') : mantissa
    return "$(mantissa)e$(exponent)"
end

# ---------------------------------------------------------------- float(), int(), round(), isclose()
const _PY_FLOAT = r"^[+-]?(?:(?:\d(?:_?\d)*)(?:\.(?:\d(?:_?\d)*)?)?|\.\d(?:_?\d)*)(?:[eE][+-]?\d(?:_?\d)*)?$"
const _PY_SPECIAL = r"^[+-]?(?:inf|infinity|nan)$"i

"""
    py_tryfloat(x) -> Float64 or nothing

Python `float(x)` without the exception: numbers (and bools) pass through,
strings follow Python's grammar (surrounding whitespace, underscores between
digits, inf/infinity/nan in any case). Anything else gives nothing.
"""
py_tryfloat(x::Real) = Float64(x)
function py_tryfloat(x::AbstractString)
    s = strip(x)
    if occursin(_PY_FLOAT, s)
        return parse(Float64, replace(s, '_' => ""))
    elseif occursin(_PY_SPECIAL, s)
        v = occursin(r"nan"i, s) ? NaN : Inf
        return startswith(s, '-') ? -v : v
    end
    return nothing
end
py_tryfloat(x) = nothing

"""
    py_float(x) -> Float64

Python `float(x)`; throws ArgumentError with Python's ValueError/TypeError text.
"""
function py_float(x)
    v = py_tryfloat(x)
    v === nothing || return v
    x isa AbstractString && throw(ArgumentError("could not convert string to float: $(py_repr(x))"))
    throw(ArgumentError("float() argument must be a string or a real number, not '$(py_typename(x))'"))
end

py_typename(::Nothing) = "NoneType"
py_typename(::Missing) = "float"   # pandas NaN
py_typename(::AbstractVector) = "list"
py_typename(::AbstractDict) = "dict"
py_typename(x) = string(typeof(x))

"""
    py_int(x) -> Int

Python `int(x)` of a string (decimal digits, optional sign, underscores) or a
number (floats truncated); ArgumentError with Python's message otherwise.
"""
function py_int(x::AbstractString)
    s = strip(x)
    occursin(r"^[+-]?\d(?:_?\d)*$", s) && return parse(Int, replace(s, '_' => ""))
    throw(ArgumentError("invalid literal for int() with base 10: $(py_repr(String(x)))"))
end
py_int(x::Integer) = Int(x)
py_int(x::AbstractFloat) = isfinite(x) ? Int(trunc(x)) :
    throw(ArgumentError(isnan(x) ? "cannot convert float NaN to integer" : "cannot convert float infinity to integer"))
py_int(::Nothing) = throw(ArgumentError("int() argument must be a string, a bytes-like object or a real number, not 'NoneType'"))
py_int(x) = throw(ArgumentError("int() argument must be a string, a bytes-like object or a real number, not '$(py_typename(x))'"))

"""
    py_round(x, n) -> Float64

Python's `round(x, n)` for a float: the decimal value correctly rounded to `n`
places (ties to even on the exact binary value), read back as a float.
"""
function py_round(x::Real, n::Integer)
    x = Float64(x)
    isfinite(x) || return x
    return parse(Float64, Printf.format(Printf.Format("%.$(n)f"), x))
end

"""Python's `math.isclose(a, b; rel_tol=1e-9, abs_tol=0.0)`."""
function py_isclose(a, b; rel_tol=1e-9, abs_tol=0.0)
    a, b = Float64(a), Float64(b)
    a == b && return true
    (isinf(a) || isinf(b)) && return false
    diff = abs(b - a)
    return diff <= abs(rel_tol * b) || diff <= abs(rel_tol * a) || diff <= abs_tol
end

"""Python truthiness of a value (None, False, 0, "", empty containers are false)."""
py_truth(::Nothing) = false
py_truth(x::Bool) = x
py_truth(x::Number) = !iszero(x)
py_truth(x::AbstractString) = !isempty(x)
py_truth(x::Union{AbstractVector,AbstractDict,Tuple,AbstractSet}) = !isempty(x)
py_truth(x) = true

"""`sorted(x)` of strings, as Python sorts them (by code point)."""
py_sorted(x) = sort!(collect(x))

# ---------------------------------------------------------------- strings and text files
"""
    py_slice(s, a, b)

Python `s[a:b]` for 0 <= a <= b on characters (not bytes), clamped to the string.
"""
function py_slice(s::AbstractString, a::Integer, b::Integer)
    if isascii(s)
        n = ncodeunits(s)
        lo, hi = min(a, n), min(b, n)
        return lo >= hi ? "" : String(codeunits(s)[lo+1:hi])
    end
    c = collect(s)
    lo, hi = min(a, length(c)), min(b, length(c))
    return lo >= hi ? "" : String(c[lo+1:hi])
end

"""
    read_text(path) -> String

The file decoded as UTF-8 dropping invalid bytes (Python `errors='ignore'`),
a leading BOM kept as Python does, and newlines translated as Python's text
mode does (`\\r\\n` and `\\r` become `\\n`).
"""
function read_text(path)
    raw = read(path, String)
    text = isvalid(raw) ? raw : String(filter(isvalid, collect(raw)))
    return replace(text, "\r\n" => "\n", '\r' => '\n')
end

"""
    py_splitlines(text; keepends=false) -> Vector{String}

Python `str.splitlines`: breaks at \\n, \\r, \\r\\n, \\v, \\f, \\x1c-\\x1e,
\\x85, \\u2028, \\u2029.
"""
function py_splitlines(text::AbstractString; keepends::Bool=false)
    out = String[]
    buf = IOBuffer()
    chars = collect(text)
    i = 1
    while i <= length(chars)
        c = chars[i]
        if c in ('\n', '\r', '\v', '\f', '\x1c', '\x1d', '\x1e', '\u85', ' ', ' ')
            ending = (c == '\r' && i < length(chars) && chars[i+1] == '\n') ? "\r\n" : string(c)
            keepends && print(buf, ending)
            push!(out, String(take!(buf)))
            i += length(ending)
        else
            print(buf, c)
            i += 1
        end
    end
    rest = String(take!(buf))
    isempty(rest) || push!(out, rest)
    return out
end
