# Small pieces of Python behaviour the data layer depends on, so that the Julia
# port reads the same files into the same values and error messages:
#
#   py_float, py_int      Python float(str) / int(str), with Python's messages
#   py_repr               repr() of a float or str (shortest round trip, Python layout)
#   casefold              str.casefold()
#   py_slice              s[a:b] on characters (0-based, clamped)
#   read_text             open(path, encoding='utf-8', errors='ignore').read() (universal newlines)
#   py_splitlines         str.splitlines(keepends)
#   parsefile_ordered     tomllib.load: TOML with the file's key order (OrderedDict)
#
# Python ValueError/TypeError become ArgumentError with the same message,
# KeyError becomes KeyError(message).
#
# These reproduce Python's *behaviour* (what is read, which texts and errors
# result). Arithmetic is plain Julia: results are compared with Python within
# tolerances, not bit for bit.

isna(x) = x === missing || x === nothing || (x isa AbstractFloat && isnan(x))   # pd.isna

casefold(s::AbstractString) = Base.Unicode.normalize(String(s); casefold=true, compose=false)

const _PY_FLOAT = r"^[+-]?(?:(?:\d(?:_?\d)*)(?:\.(?:\d(?:_?\d)*)?)?|\.\d(?:_?\d)*)(?:[eE][+-]?\d(?:_?\d)*)?$"
const _PY_SPECIAL = r"^[+-]?(?:inf|infinity|nan)$"i

"""
    py_tryfloat(x) -> Float64 or nothing

Python `float(x)` without the exception: numbers pass through, strings follow
Python's grammar (surrounding whitespace, underscores between digits,
inf/infinity/nan in any case). Anything else (missing, nothing) gives nothing.
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

Python `float(x)`; raises ArgumentError with Python's ValueError/TypeError text.
"""
function py_float(x)
    v = py_tryfloat(x)
    v === nothing || return v
    x isa AbstractString && throw(ArgumentError("could not convert string to float: $(py_repr(x))"))
    throw(ArgumentError("float() argument must be a string or a real number, not '$(py_typename(x))'"))
end

py_typename(::Nothing) = "NoneType"
py_typename(::Missing) = "float"   # pandas NaN
py_typename(x) = string(typeof(x))

"""
    py_int(x) -> Int

Python `int(x)` of a string (decimal digits, optional sign, underscores) or an
integral number; ArgumentError with Python's message otherwise.
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

"""
    py_repr(x)

Python `repr` of a float (shortest round-trip digits, Python's layout: `1e-05`,
`1e+16`, `100.0`) or a str (single quotes unless the text contains one and no
double quote). `py_str` is `str()`: the same for floats, the text itself for strings.
"""
function py_repr(x::AbstractFloat)
    x = Float64(x)
    isnan(x) && return "nan"
    isinf(x) && return x > 0 ? "inf" : "-inf"
    x == 0 && return signbit(x) ? "-0.0" : "0.0"
    s = string(x)                      # Julia: shortest round trip, e.g. "1.0e-5"
    sign = startswith(s, '-') ? "-" : ""
    s = lstrip(s, '-')
    mant, ex = occursin('e', s) ? split(s, 'e') : (s, "0")
    ip, fp = occursin('.', mant) ? split(mant, '.') : (mant, "")
    digits = ip * fp
    point = length(ip) + parse(Int, ex)          # position of the decimal point
    lead = length(digits) - length(lstrip(digits, '0'))
    digits = rstrip(lstrip(digits, '0'), '0')
    point -= lead
    n = length(digits)
    exp10 = point - 1                            # scientific exponent
    if -4 <= exp10 < 16
        if point <= 0
            body = "0." * "0"^(-point) * digits
        elseif point >= n
            body = digits * "0"^(point - n) * ".0"
        else
            body = digits[1:point] * "." * digits[point+1:end]
        end
    else
        m = n == 1 ? digits : digits[1:1] * "." * digits[2:end]
        body = m * "e" * (exp10 < 0 ? "-" : "+") * lpad(string(abs(exp10)), 2, '0')
    end
    return sign * body
end
function py_repr(s::AbstractString)
    q = occursin('\'', s) && !occursin('"', s) ? '"' : '\''
    body = replace(String(s), '\\' => "\\\\", '\n' => "\\n", '\r' => "\\r", '\t' => "\\t")
    q == '\'' && (body = replace(body, '\'' => "\\'"))
    return string(q, body, q)
end
py_repr(::Nothing) = "None"
py_repr(x::Integer) = string(x)
py_repr(x::Bool) = x ? "True" : "False"

py_str(x::AbstractString) = String(x)
py_str(x::AbstractFloat) = py_repr(x)
py_str(x::Integer) = string(x)
py_str(::Nothing) = "None"
py_str(::Missing) = "nan"
py_str(x) = string(x)

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

# ---------------------------------------------------------------- TOML in file order
# Julia's TOML stdlib returns Dicts, which lose the key order; Python's tomllib
# keeps it, and the registry depends on it (species order). The values come
# from the stdlib parser; this scanner only records where each key first
# appears, and the tables are rebuilt as OrderedDicts in that order.

"""
    parsefile_ordered(path) / parse_ordered(text)

TOML like `TOML.parsefile`, with every table an `OrderedDict{String,Any}` in
the order of the file (as Python's `tomllib`).
"""
parsefile_ordered(path::AbstractString) = parse_ordered(read(path, String))
function parse_ordered(text::AbstractString)
    data = TOML.parse(text)
    order = _toml_key_order(text)
    return _ordered(data, order, String[])
end

function _ordered(d::AbstractDict, order, path)
    keys_ = collect(keys(d))
    position(k) = get(order, vcat(path, k), typemax(Int))
    sort!(keys_; by=k -> (position(k), k))
    out = OrderedDict{String,Any}()
    for k in keys_
        out[k] = _ordered(d[k], order, vcat(path, k))
    end
    return out
end
_ordered(v::AbstractVector, order, path) = Any[_ordered(x, order, path) for x in v]
_ordered(v, order, path) = v

function _toml_key_order(text::AbstractString)
    s = collect(text)
    n = length(s)
    order = Dict{Vector{String},Int}()
    counter = Ref(0)
    note(path) = for k in 1:length(path)
        p = path[1:k]
        haskey(order, p) || (order[p] = (counter[] += 1))
    end
    i = Ref(1)
    peek(k=0) = i[] + k <= n ? s[i[]+k] : '\0'
    function skip_ws(newlines::Bool)
        while i[] <= n
            c = s[i[]]
            if c == ' ' || c == '\t' || (newlines && (c == '\n' || c == '\r'))
                i[] += 1
            elseif c == '#'
                while i[] <= n && s[i[]] != '\n'
                    i[] += 1
                end
            else
                break
            end
        end
    end
    function read_string()           # at a quote; returns the (unescaped enough) text
        q = s[i[]]
        if peek(1) == q && peek(2) == q   # multi-line
            i[] += 3
            start = i[]
            while i[] <= n && !(s[i[]] == q && peek(1) == q && peek(2) == q)
                q == '"' && s[i[]] == '\\' && (i[] += 1)
                i[] += 1
            end
            str = String(s[start:min(i[]-1, n)])
            i[] += 3
            while peek() == q   # up to two extra quotes belong to the string
                i[] += 1
            end
            return str
        end
        i[] += 1
        buf = IOBuffer()
        while i[] <= n && s[i[]] != q
            if q == '"' && s[i[]] == '\\'
                i[] += 1
                e = peek()
                if e == 'u' || e == 'U'
                    len = e == 'u' ? 4 : 8
                    print(buf, Char(parse(UInt32, String(s[i[]+1:i[]+len]); base=16)))
                    i[] += len
                else
                    print(buf, get(Dict('n' => '\n', 't' => '\t', 'r' => '\r', 'b' => '\b',
                                        'f' => '\f', '"' => '"', '\\' => '\\', 'e' => '\e'), e, e))
                end
            else
                print(buf, s[i[]])
            end
            i[] += 1
        end
        i[] += 1
        return String(take!(buf))
    end
    function read_key()
        parts = String[]
        while true
            skip_ws(false)
            c = peek()
            if c == '"' || c == '\''
                push!(parts, read_string())
            else
                start = i[]
                while i[] <= n && (isletter(s[i[]]) || isdigit(s[i[]]) || s[i[]] in ('_', '-'))
                    i[] += 1
                end
                push!(parts, String(s[start:i[]-1]))
            end
            skip_ws(false)
            peek() == '.' || break
            i[] += 1
        end
        return parts
    end
    function skip_value(path)
        skip_ws(false)
        c = peek()
        if c == '"' || c == '\''
            read_string()
        elseif c == '['
            i[] += 1
            while true
                skip_ws(true)
                peek() == ']' && (i[] += 1; break)
                i[] > n && break
                skip_value(path)
                skip_ws(true)
                peek() == ',' && (i[] += 1)
            end
        elseif c == '{'
            i[] += 1
            while true
                skip_ws(false)
                peek() == '}' && (i[] += 1; break)
                i[] > n && break
                key = read_key()
                note(vcat(path, key))
                skip_ws(false)
                peek() == '=' && (i[] += 1)
                skip_value(vcat(path, key))
                skip_ws(false)
                peek() == ',' && (i[] += 1)
            end
        else
            while i[] <= n && !(s[i[]] in (',', ']', '}', '\n', '\r', '#'))
                i[] += 1
            end
        end
    end
    table = String[]
    while true
        skip_ws(true)
        i[] > n && break
        if peek() == '['
            double = peek(1) == '['
            i[] += double ? 2 : 1
            table = read_key()
            note(table)
            skip_ws(false)
            i[] += double ? 2 : 1
        else
            key = read_key()
            if all(isempty, key) && peek() != '='   # not a key (cannot happen in valid TOML)
                i[] += 1
                continue
            end
            note(vcat(table, key))
            skip_ws(false)
            peek() == '=' && (i[] += 1)
            skip_value(vcat(table, key))
        end
    end
    return order
end
