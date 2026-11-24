# Small helpers that reproduce Python behaviour the ported modules rely on:
# `ValueError` (the Python exception family most errors belong to), `pyrepr`
# (Python's repr, used inside error messages), `pyfloat` (Python's float() of a
# number or a string), `pyround` (Python's round(x, n)) and `pyisclose`
# (math.isclose).

"""
Python's `ValueError` family. Error types of the ported modules that are
ValueErrors in Python (`CaseError`, `FeedError`, `LegacyConversionError`)
subtype it, so `e isa ValueError` catches what `except ValueError` catches.
"""
abstract type ValueError <: Exception end
Base.showerror(io::IO, e::ValueError) = print(io, e.msg)

# ---------------------------------------------------------------- repr
"""
    pyrepr(x) -> String

Python's `repr` for the values case files hold: None, bools, ints, floats
(shortest round-trip digits, Python's switch to exponent notation), strings
(Python's quoting), lists/vectors, tuples and dicts.
"""
pyrepr(::Nothing) = "None"
pyrepr(x::Bool) = x ? "True" : "False"
pyrepr(x::Integer) = string(x)
pyrepr(x::AbstractFloat) = _pyfloat_repr(Float64(x))
pyrepr(x::AbstractString) = _pystr_repr(String(x))
pyrepr(x::Symbol) = _pystr_repr(String(x))
pyrepr(x::AbstractVector) = "[" * join((pyrepr(v) for v in x), ", ") * "]"
pyrepr(x::Tuple) = length(x) == 1 ? "(" * pyrepr(x[1]) * ",)" : "(" * join((pyrepr(v) for v in x), ", ") * ")"
pyrepr(x::AbstractDict) = "{" * join((pyrepr(k) * ": " * pyrepr(v) for (k, v) in x), ", ") * "}"
pyrepr(x) = repr(x)

"""Python's `str()`: like `pyrepr`, except that strings are not quoted."""
pystr(x::AbstractString) = String(x)
pystr(x) = pyrepr(x)

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
    pyformat_g(value, p) -> String

Python's `f"{value:.{p}g}"` (C's "%.<p>g"). Julia's Printf `%g` picks fixed or
exponent notation before rounding (999999.5 with p = 6 gives "1000000"; C and
Python give "1e+06"), so the exponent is taken from the rounded "%e" form, as the
C standard specifies. Trailing zeros are removed.
"""
function pyformat_g(value::Real, p::Integer)
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

# ---------------------------------------------------------------- float(), round(), isclose()
const _PYFLOAT_TEXT = r"^[+-]?(?:(?:\d(?:_?\d)*)?\.?(?:\d(?:_?\d)*)?(?:[eE][+-]?\d(?:_?\d)*)?|inf|infinity|nan)$"i

"""
    pyfloat(x) -> Float64

Python's `float(x)` for numbers, bools and strings ("1_000", " 2.5 ", "inf").
Anything else, or a string Python would refuse, throws an `ArgumentError`
(Python: TypeError or ValueError); callers turn it into their own message.
"""
pyfloat(x::Real) = Float64(x)
function pyfloat(x::AbstractString)
    s = strip(x)
    if !isempty(s) && occursin(_PYFLOAT_TEXT, s) && any(c -> isdigit(c) || c in "iInN", s)
        v = tryparse(Float64, replace(s, '_' => ""))
        v === nothing || return v
    end
    throw(ArgumentError("could not convert string to float: $(pyrepr(String(x)))"))
end
pyfloat(x) = throw(ArgumentError("float() argument must be a string or a real number, not $(typeof(x))"))

"""
    pyround(x, n) -> Float64

Python's `round(x, n)` for a float: the decimal value correctly rounded to `n`
places (ties to even on the exact binary value), read back as a float.
"""
function pyround(x::Real, n::Integer)
    x = Float64(x)
    isfinite(x) || return x
    return parse(Float64, Printf.format(Printf.Format("%.$(n)f"), x))
end

"""Python's `math.isclose(a, b; rel_tol=1e-9, abs_tol=0.0)`."""
function pyisclose(a, b; rel_tol=1e-9, abs_tol=0.0)
    a, b = Float64(a), Float64(b)
    a == b && return true
    (isinf(a) || isinf(b)) && return false
    diff = abs(b - a)
    return diff <= abs(rel_tol * b) || diff <= abs(rel_tol * a) || diff <= abs_tol
end

"""Python truthiness of a value (None, False, 0, "", empty containers are false)."""
pytruth(::Nothing) = false
pytruth(x::Bool) = x
pytruth(x::Number) = !iszero(x)
pytruth(x::AbstractString) = !isempty(x)
pytruth(x::Union{AbstractVector,AbstractDict,Tuple,AbstractSet}) = !isempty(x)
pytruth(x) = true

"""`sorted(x)` of strings, as Python sorts them (by code point)."""
pysorted(x) = sort!(collect(x))
