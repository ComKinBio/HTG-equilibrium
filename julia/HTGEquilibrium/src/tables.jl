# Small pandas-like tables for the result tables (Python: the pandas DataFrames
# of htg/results.py, `DataFrame.to_csv(index=False)` and `pd.read_csv`).
#
# The type is the data layer's `DataTable` (databases.jl). For result tables a
# column holds one Julia type per pandas dtype:
#
#     float64   Float64 (NaN where pandas has NaN)       numeric
#     int64     Int (only when no value is missing)      numeric
#     bool      Bool (only when no value is missing)
#     object    String, Bool, Int, Float64 or `missing` (also a column of only None)
#
# `DataTable(rows)` infers these from a list of dicts as `pd.DataFrame(rows)`
# does (columns in order of first appearance; an int column with a missing
# value becomes float64, a bool column with a missing value object), so that
# `to_csv` writes what pandas writes: floats as Python's repr ("300.0",
# "1e-06"), NaN and missing as "", bools as True/False, text quoted only when
# it holds a comma, a quote or a line break (csv.QUOTE_MINIMAL). `read_csv`
# reads such files back with pandas' type inference.

_isnull(x) = x === missing || x === nothing || (x isa AbstractFloat && isnan(x))

"""
Values of one column as pandas stores them: (values, numeric). Ints with a
missing value or mixed with floats become Float64 (NaN for missing).
"""
function _infer_column(values::AbstractVector)
    nonnull = [v for v in values if !_isnull(v)]
    has_null = length(nonnull) < length(values)
    if isempty(nonnull) && !isempty(values) && all(v -> v === nothing, values)
        return Any[missing for _ in values], false   # only None: object (pandas)
    elseif all(v -> v isa Real && !(v isa Bool), nonnull)
        if !has_null && !isempty(nonnull) && all(v -> v isa Integer, nonnull)
            return Any[Int(v) for v in values], true
        end
        return Any[_isnull(v) ? NaN : Float64(v) for v in values], true
    elseif !has_null && all(v -> v isa Bool, nonnull)
        return Any[v for v in values], false
    end
    return Any[_isnull(v) ? missing : (v isa AbstractString ? String(v) : v) for v in values], false
end

function _table(columns, data)
    out = DataTable(String[], Dict{String,Vector{Any}}(), Set{String}(), Int[])
    n = isempty(columns) ? 0 : length(data[first(columns)])
    for c in columns
        values, numeric = _infer_column(data[c])
        push!(out.columns, c)
        out.data[c] = values
        numeric && push!(out.numeric, c)
    end
    out.index = collect(0:n-1)
    return out
end

"""
    DataTable(rows; columns=nothing)

`pd.DataFrame(rows)` of a list of dicts: the columns in order of first
appearance (or `columns`), a missing key is NaN; column types as pandas infers them.
"""
function DataTable(rows::AbstractVector{<:AbstractDict}; columns=nothing)
    names = columns === nothing ? String[] : collect(String, columns)
    if columns === nothing
        seen = Set{String}()
        for row in rows, k in keys(row)
            k = String(k)
            k in seen || (push!(seen, k); push!(names, k))
        end
    end
    data = Dict{String,Vector{Any}}(c => Any[get(row, c, missing) for row in rows] for c in names)
    return _table(names, data)
end

"""An empty table (no columns, no rows): `pd.DataFrame()`."""
DataTable() = DataTable(OrderedDict{String,Any}[])

Base.getindex(t::DataTable, c::AbstractString) =
    haskey(t.data, c) ? t.data[c] : throw(KeyError(c))
Base.haskey(t::DataTable, c::AbstractString) = haskey(t.data, c)
"""pandas `DataFrame.empty`: no rows or no columns."""
Base.isempty(t::DataTable) = isempty(t.columns) || nrows(t) == 0
Base.copy(t::DataTable) = DataTable(copy(t.columns), Dict{String,Vector{Any}}(c => copy(v) for (c, v) in t.data),
                                    copy(t.numeric), copy(t.index))
Base.size(t::DataTable) = (nrows(t), length(t.columns))

"""Row `i` (1-based) as an ordered dict of column => value."""
row(t::DataTable, i::Integer) = OrderedDict{String,Any}(c => t.data[c][i] for c in t.columns)
"""All rows as ordered dicts (pandas `to_dict('records')`)."""
rows(t::DataTable) = [row(t, i) for i in 1:nrows(t)]

"""
Set column `c` to `values` (a new column goes last; pandas `t[c] = values`).
A single value is repeated for every row. Types are inferred again.
"""
function setcolumn!(t::DataTable, c::AbstractString, values)
    c = String(c)
    n = isempty(t.columns) ? (values isa AbstractVector ? length(values) : 0) : nrows(t)
    v = values isa AbstractVector ? Any[x for x in values] : Any[values for _ in 1:n]
    length(v) == n || throw(ArgumentError("Length of values ($(length(v))) does not match length of index ($n)"))
    col, numeric = _infer_column(v)
    c in t.columns || push!(t.columns, c)
    t.data[c] = col
    numeric ? push!(t.numeric, c) : delete!(t.numeric, c)
    isempty(t.index) && n > 0 && (t.index = collect(0:n-1))
    return t
end

"""pandas `t.insert(position, c, values)`: position is 0-based, as in pandas."""
function insert_column!(t::DataTable, position::Integer, c::AbstractString, values)
    c in t.columns && throw(ArgumentError("cannot insert $c, already exists"))
    setcolumn!(t, c, values)
    pop!(t.columns)
    insert!(t.columns, position + 1, String(c))
    return t
end

"""The table with only `columns`, in that order."""
select(t::DataTable, columns) = DataTable(collect(String, columns),
                                          Dict{String,Vector{Any}}(c => copy(t[c]) for c in columns),
                                          Set{String}(c for c in columns if c in t.numeric), copy(t.index))

"""The rows at positions `rows` (1-based; a Bool mask also works), keeping the index labels."""
function take_rows(t::DataTable, rows::AbstractVector)
    idx = eltype(rows) == Bool ? findall(rows) : collect(Int, rows)
    return DataTable(copy(t.columns), Dict{String,Vector{Any}}(c => t.data[c][idx] for c in t.columns),
                     copy(t.numeric), t.index[idx])
end

"""
    concat(tables) -> DataTable

`pd.concat(tables, ignore_index=True)`: the rows of all tables, the columns of
all tables in order of first appearance; column types inferred again.
"""
function concat(tables)
    names = String[]
    for t in tables, c in t.columns
        c in names || push!(names, c)
    end
    data = Dict{String,Vector{Any}}(c => Any[] for c in names)
    for t in tables, c in names
        append!(data[c], haskey(t.data, c) ? t.data[c] : Any[missing for _ in 1:nrows(t)])
    end
    return _table(names, data)
end

"""Equal columns and values; NaN equals NaN and missing equals missing."""
Base.isequal(a::DataTable, b::DataTable) =
    a.columns == b.columns && all(isequal(a.data[c], b.data[c]) for c in a.columns)

function Base.show(io::IO, t::DataTable)
    print(io, "DataTable($(nrows(t)) rows: ", join(t.columns, ", "), ")")
end

# ---------------------------------------------------------------- CSV (pandas to_csv / read_csv)
_csv_cell(::Union{Missing,Nothing}) = ""
_csv_cell(x::Bool) = x ? "True" : "False"
_csv_cell(x::Integer) = string(x)
_csv_cell(x::AbstractFloat) = isnan(x) ? "" : py_repr(x)
_csv_cell(x::AbstractString) = String(x)
_csv_cell(x) = py_str(x)

_csv_quote(s::String) =
    any(c -> c == ',' || c == '"' || c == '\n' || c == '\r', s) ? "\"" * replace(s, "\"" => "\"\"") * "\"" : s

function _csv_line(io, fields)
    # Python's csv writer quotes a row that is a single empty field
    if length(fields) == 1 && isempty(fields[1])
        print(io, "\"\"\n")
    else
        join(io, (_csv_quote(f) for f in fields), ',')
        print(io, '\n')
    end
end

"""The text `DataFrame.to_csv(index=False)` writes (pandas 3, Python's csv module)."""
function csv_text(t::DataTable)
    io = IOBuffer()
    if isempty(t.columns)
        print(io, '\n')
        return String(take!(io))
    end
    _csv_line(io, t.columns)
    for i in 1:nrows(t)
        _csv_line(io, [_csv_cell(t.data[c][i]) for c in t.columns])
    end
    return String(take!(io))
end

"""Write `t` as `DataFrame.to_csv(path, index=False)` does."""
to_csv(t::DataTable, path::AbstractString) = (write(path, csv_text(t)); nothing)

const _PANDAS_INT = r"^[+-]?\d+$"
const _PANDAS_TRUE = ("True", "TRUE", "true")
const _PANDAS_FALSE = ("False", "FALSE", "false")

"""
    read_csv(path) -> DataTable

A results CSV as `pd.read_csv` reads it (Python: `htg.plotting.read_csv`): the
default NA strings are missing; a column of integers without missing values is
Int, a column of numbers Float64 (NaN where empty), a column of True/False
without missing values Bool, anything else text (`missing` where empty). An
empty file gives an empty table.
"""
function read_csv(path::AbstractString)
    isempty(strip(read(path, String))) && return DataTable()
    raw, n = _csv_strings(path)
    out = DataTable(String[], Dict{String,Vector{Any}}(), Set{String}(), collect(0:n-1))
    for (name, values) in raw
        values = Any[(v === missing || v in PANDAS_NA) ? missing : v for v in values]
        nonnull = [v for v in values if v !== missing]
        has_null = length(nonnull) < n
        if !has_null && !isempty(nonnull) && all(v -> occursin(_PANDAS_INT, v), nonnull)
            col = Any[parse(Int, v) for v in values]
            push!(out.numeric, name)
        elseif all(v -> occursin(PANDAS_NUMBER, v), nonnull)
            col = Any[v === missing ? NaN : py_float(strip(v)) for v in values]
            push!(out.numeric, name)
        elseif all(v -> v in _PANDAS_TRUE || v in _PANDAS_FALSE, nonnull)
            col = Any[v === missing ? missing : v in _PANDAS_TRUE for v in values]
        else
            col = values
        end
        push!(out.columns, name)
        out.data[name] = col
    end
    return out
end

# ---------------------------------------------------------------- text (pandas to_string)
# `DataFrame.to_string(index=False)` of pandas 3.1.1 with the display options of
# the command line (`htg.cli._print_table`: no row/column limits; to_string
# itself never truncates a value). pandas formats each column by its dtype
# (format.py: FloatArrayFormatter, _IntArrayFormatter, _GenericArrayFormatter,
# _ExtensionArrayFormatter for the `str` dtype), right-justifies values and
# header to the widest of them, puts a space before the header of a numeric
# column (int, float and bool), and joins the columns with one space. Widths
# count characters (Python's len), not display width.

"""The pandas dtype a table column would have: :int, :float, :bool, :str or :object."""
function _pd_kind(t::DataTable, c)
    values = t[c]
    c in t.numeric && return all(v -> v isa Integer && !(v isa Bool), values) ? :int : :float
    !isempty(values) && all(v -> v isa Bool, values) && return :bool
    nonnull = [v for v in values if !_isnull(v)]
    !isempty(nonnull) && all(v -> v isa AbstractString, nonnull) && return :str
    return :object
end

"""pprint_thing with escape_chars ("\\t", "\\r", "\\n")."""
_pd_escape(s) = replace(String(s), "\t" => "\\t", "\r" => "\\r", "\n" => "\\n")

"""Python's `format(x, '.6f')` / `format(x, '.6e')` (inf and nan in lower case)."""
function _pd_fmt(x::Float64, style::Char)
    isnan(x) && return "nan"
    isinf(x) && return x > 0 ? "inf" : "-inf"
    return style == 'f' ? @sprintf("%.6f", x) : @sprintf("%.6e", x)
end

const _PD_DECIMAL = r"^\s*[+-]?[0-9]+\.[0-9]*$"

"""pandas `_trim_zeros_float`: drop trailing zeros equally from all decimal numbers, keep one."""
function _pd_trim_zeros(strings::Vector{String})
    isnum = [occursin(_PD_DECIMAL, s) for s in strings]
    any(isnum) || return strings
    out = copy(strings)
    while all(endswith(s, "0") for (s, n) in zip(out, isnum) if n)
        out = [n ? String(chop(s)) : s for (s, n) in zip(out, isnum)]
    end
    return [n && endswith(s, ".") ? s * "0" : s for (s, n) in zip(out, isnum)]
end

"""pandas `FloatArrayFormatter` (precision 6, no leading space): fixed, or scientific when needed."""
function _pd_float_column(values)
    x = Float64[Float64(v) for v in values]
    formatted(style) = _pd_trim_zeros(String[isnan(v) ? "NaN" : _pd_fmt(v, style) for v in x])
    fixed = formatted('f')
    too_long = !isempty(fixed) && maximum(length, fixed) > 6 + 6
    large = any(v -> abs(v) > 1e6, x)
    small = any(v -> 0 < abs(v) < 1e-6, x)
    return (small || (too_long && large)) ? formatted('e') : fixed
end

"""A value of an object column (`_GenericArrayFormatter`): None, NaN; floats as `f"{x: .6f}"` trimmed."""
function _pd_object_cell(v)
    v === nothing && return "None"
    _isnull(v) && return "NaN"
    v isa Bool && return v ? "True" : "False"
    if v isa AbstractFloat
        s = rstrip((signbit(v) ? "" : " ") * _pd_fmt(Float64(v), 'f'), '0')
        return endswith(s, ".") ? s * "0" : String(s)
    end
    return _pd_escape(v isa AbstractString ? v : py_str(v))
end

function _pd_column(t::DataTable, c)
    kind = _pd_kind(t, c)
    values = t[c]
    cells = if kind == :float
        _pd_float_column(values)
    elseif kind == :int
        String[string(v) for v in values]
    elseif kind == :bool
        String[v ? "True" : "False" for v in values]
    elseif kind == :str
        String[_isnull(v) ? "NaN" : _pd_escape(v) for v in values]
    else
        String[_pd_object_cell(v) for v in values]
    end
    header = kind in (:float, :int, :bool) ? " " * c : c
    return header, cells
end

_rjust(s, n) = " "^max(0, n - length(s)) * s

"""
    to_string(table) -> String

The text of pandas `DataFrame.to_string(index=False)`, as the command line
prints tables: values and headers right-justified, numeric headers (int, float,
bool) with a leading space, floats in pandas' common format per column
(`300.0`, `-481234.56789`, `1.000000e-07`), `NaN` for missing numbers and text,
`None` for Python's None (`nothing`) in an object column, bools as True/False.
A table without rows gives pandas' "Empty DataFrame" text.
"""
function to_string(t::DataTable)
    if isempty(t)
        return "Empty DataFrame\nColumns: [$(join(t.columns, ", "))]\nIndex: []"
    end
    columns = Vector{String}[]
    for c in t.columns
        header, cells = _pd_column(t, c)
        width = max(length(header), maximum(length, cells))
        push!(columns, [_rjust(s, width) for s in [header; cells]])
    end
    return join((join((col[i] for col in columns), ' ') for i in 1:nrows(t)+1), '\n')
end
