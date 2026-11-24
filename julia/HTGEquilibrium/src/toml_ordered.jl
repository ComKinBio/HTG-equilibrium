# TOML with the key order of the file (Python's tomllib keeps it; Julia's TOML
# standard library returns plain Dicts). The values come from `TOML.parse`, so
# the TOML grammar, value types and syntax errors are the standard library's;
# a light scanner of the text only records in which order the keys appear in
# each table, and the result is rebuilt as nested OrderedDicts in that order.
#
#     parse_toml(text)  -> OrderedDict{String,Any}   (Python: tomllib.loads)
#     read_toml(path)   -> OrderedDict{String,Any}   (Python: tomllib.load)

"""
    parse_toml(text) -> OrderedDict{String,Any}

TOML text as nested `OrderedDict{String,Any}`s whose keys are in the order of
the file, like Python's `tomllib.loads`. Integers stay `Int64`, floats
`Float64`. Throws `TOML.ParserError` for invalid TOML.
"""
function parse_toml(text::AbstractString)
    data = TOML.parse(String(text))
    order = _toml_key_order(String(text))
    return _ordered(data, Any[], order)
end

"""`parse_toml` of a file (read as UTF-8)."""
read_toml(path::AbstractString) = parse_toml(read(path, String))

function _ordered(x::AbstractDict, path::Vector{Any}, order)
    seen = get(order, path, String[])
    keys_in_order = [k for k in seen if haskey(x, k)]
    # Keys the scanner missed (it should miss none) follow in sorted order.
    append!(keys_in_order, sort!([k for k in keys(x) if !(k in seen)]))
    out = OrderedDict{String,Any}()
    for k in keys_in_order
        out[k] = _ordered(x[k], Any[path..., k], order)
    end
    return out
end
_ordered(x::AbstractVector, path::Vector{Any}, order) =
    Any[_ordered(v, Any[path..., i], order) for (i, v) in enumerate(x)]
_ordered(x, path, order) = x

# ---------------------------------------------------------------- key-order scanner
mutable struct _Scanner
    s::String
    i::Int                                   # byte index of the next character
    order::Dict{Vector{Any},Vector{String}}  # table path -> keys in order of appearance
    aot::Dict{Vector{Any},Int}               # array-of-tables path -> elements so far
end

_peek(sc::_Scanner, k=0) = (j = sc.i; for _ in 1:k; j = nextind(sc.s, j); end; j <= ncodeunits(sc.s) ? sc.s[j] : '\0')
_next!(sc::_Scanner) = (sc.i = nextind(sc.s, sc.i))
_atend(sc::_Scanner) = sc.i > ncodeunits(sc.s)

function _record!(sc::_Scanner, path::Vector{Any}, key::String)
    keys = get!(sc.order, copy(path), String[])
    key in keys || push!(keys, key)
end

function _skip_ws!(sc::_Scanner; newlines=false, comments=false)
    while !_atend(sc)
        c = _peek(sc)
        if c == ' ' || c == '\t' || (newlines && (c == '\n' || c == '\r'))
            _next!(sc)
        elseif comments && c == '#'
            while !_atend(sc) && _peek(sc) != '\n'
                _next!(sc)
            end
        else
            break
        end
    end
end

# A basic ("...", with escapes) or literal ('...') string, single or multi-line.
function _read_string!(sc::_Scanner)
    q = _peek(sc)
    multi = _peek(sc, 1) == q && _peek(sc, 2) == q
    start = sc.i
    if multi
        for _ in 1:3; _next!(sc); end
        while !_atend(sc)
            if q == '"' && _peek(sc) == '\\'
                _next!(sc); _next!(sc)
            elseif _peek(sc) == q && _peek(sc, 1) == q && _peek(sc, 2) == q
                for _ in 1:3; _next!(sc); end
                while _peek(sc) == q   # up to two quotes may end the content
                    _next!(sc)
                end
                break
            else
                _next!(sc)
            end
        end
    else
        _next!(sc)
        while !_atend(sc) && _peek(sc) != q
            q == '"' && _peek(sc) == '\\' && _next!(sc)
            _next!(sc)
        end
        _next!(sc)
    end
    return sc.s[start:prevind(sc.s, sc.i)]
end

# One key part: bare or quoted; returns its value.
function _read_key_part!(sc::_Scanner)
    _skip_ws!(sc)
    c = _peek(sc)
    if c == '"' || c == '\''
        raw = _read_string!(sc)
        # Let the TOML parser decode the quoted key (escapes included).
        return only(keys(TOML.parse(raw * " = 0")))
    end
    start = sc.i
    while !_atend(sc) && (isletter(_peek(sc)) || isdigit(_peek(sc)) || _peek(sc) in ('_', '-'))
        _next!(sc)
    end
    return sc.s[start:prevind(sc.s, sc.i)]
end

function _read_dotted_key!(sc::_Scanner)
    parts = String[_read_key_part!(sc)]
    _skip_ws!(sc)
    while _peek(sc) == '.'
        _next!(sc)
        push!(parts, _read_key_part!(sc))
        _skip_ws!(sc)
    end
    return parts
end

# Path of a table header's keys, with the latest element of each array of tables.
function _resolve!(sc::_Scanner, base::Vector{Any}, parts::Vector{String})
    path = copy(base)
    for k in parts
        _record!(sc, path, k)
        push!(path, k)
        haskey(sc.aot, path) && push!(path, sc.aot[path])
    end
    return path
end

# key = value inside the table at `table`.
function _key_value!(sc::_Scanner, table::Vector{Any})
    parts = _read_dotted_key!(sc)
    path = _resolve!(sc, table, parts)
    _skip_ws!(sc)
    _peek(sc) == '=' && _next!(sc)
    _skip_ws!(sc)
    _value!(sc, path)
end

function _value!(sc::_Scanner, path::Vector{Any})
    c = _peek(sc)
    if c == '"' || c == '\''
        _read_string!(sc)
    elseif c == '{'
        _next!(sc)
        _skip_ws!(sc)
        while !_atend(sc) && _peek(sc) != '}'
            _key_value!(sc, path)
            _skip_ws!(sc)
            _peek(sc) == ',' && _next!(sc)
            _skip_ws!(sc)
        end
        _next!(sc)
    elseif c == '['
        _next!(sc)
        n = 0
        while true
            _skip_ws!(sc; newlines=true, comments=true)
            (_atend(sc) || _peek(sc) == ']') && break
            n += 1
            _value!(sc, Any[path..., n])
            _skip_ws!(sc; newlines=true, comments=true)
            _peek(sc) == ',' && _next!(sc)
        end
        _next!(sc)
    else  # number, bool, date: up to the next delimiter
        while !_atend(sc) && !(_peek(sc) in (',', ']', '}', '\n', '\r', '#'))
            _next!(sc)
        end
    end
end

function _toml_key_order(s::String)
    sc = _Scanner(s, 1, Dict{Vector{Any},Vector{String}}(), Dict{Vector{Any},Int}())
    table = Any[]
    while true
        _skip_ws!(sc; newlines=true, comments=true)
        _atend(sc) && break
        if _peek(sc) == '['
            array = _peek(sc, 1) == '['
            _next!(sc); array && _next!(sc)
            parts = _read_dotted_key!(sc)
            if array
                parent = _resolve!(sc, Any[], parts[1:end-1])
                _record!(sc, parent, parts[end])
                key = Any[parent..., parts[end]]
                sc.aot[key] = get(sc.aot, key, 0) + 1
                table = Any[key..., sc.aot[key]]
            else
                table = _resolve!(sc, Any[], parts)
            end
            while !_atend(sc) && _peek(sc) != '\n'  # the closing bracket(s), a comment
                _next!(sc)
            end
        else
            _key_value!(sc, table)
            _skip_ws!(sc; comments=true)
        end
    end
    return sc.order
end
