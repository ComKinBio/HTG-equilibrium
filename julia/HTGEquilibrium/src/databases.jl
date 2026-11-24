# Readers for the four thermodynamic databases (Python: htg/databases.py).
#
#     nasa      NASA-7 polynomial records (Thermdat_WasteMech.txt)
#     pcsaft    PC-SAFT pure-component parameters (PROJECT_PCSAFT.csv)
#     epcsaft   ePC-SAFT species parameters, neutral and ionic (PROJECT_EPCSAFT.csv)
#     supcrt    SUPCRT92 aqueous HKF and mineral records (supcrt92.xml)
#
# The readers only parse files; interpreting the records is done in
# standard_states.jl, saft.jl and registry.jl. Everything here is built once by
# `load_databases` (and the user layer, before the model is shared) and read
# only afterwards, so the structures can be shared between threads.

const DATABASE_FILES = OrderedDict(
    "nasa" => "Thermdat_WasteMech.txt",
    "pcsaft" => "PROJECT_PCSAFT.csv",
    "epcsaft" => "PROJECT_EPCSAFT.csv",
    "supcrt" => "supcrt92.xml",
)

# ---------------------------------------------------------------- CSV tables
"""
    DataTable

A SAFT CSV table as pandas reads it. `columns` in file order; `data[column]`
holds one value per row: in `numeric` columns (pandas float64/int64) a Float64,
NaN when empty; in the other columns a String or `missing` (pandas NaN).
`index` holds the pandas index labels: the row number in the file (0-based),
negative for rows of the user layer, which come first.
"""
mutable struct DataTable
    columns::Vector{String}
    data::Dict{String,Vector{Any}}
    numeric::Set{String}
    index::Vector{Int}
end

nrows(t::DataTable) = length(t.index)
hascolumn(t::DataTable, c::AbstractString) = haskey(t.data, c)
"""Value in column `c` at row position `i` (1-based); `missing` if there is no such column (Python `row.get`)."""
cell(t::DataTable, c::AbstractString, i::Integer) = haskey(t.data, c) ? t.data[c][i] : missing
"""Rows (positions) where column `c` as text equals `value` (pandas `astype(str).eq`; NaN never matches)."""
rows_where(t::DataTable, c, value; rows=1:nrows(t)) =
    [i for i in rows if !isna(cell(t, c, i)) && py_str(cell(t, c, i)) == value]

# pandas' default NA strings (read_csv na_values)
const PANDAS_NA = Set(["", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND",
                       "1.#QNAN", "<NA>", "N/A", "NA", "NULL", "NaN", "None", "n/a", "nan", "null"])
const PANDAS_NUMBER = r"^\s*[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?\s*$|^\s*[+-]?(?:inf|infinity)\s*$"i

function _csv_strings(path)
    file = CSV.File(path; types=String, normalizenames=false)
    names = [lstrip(String(n), '﻿') for n in propertynames(file)]
    columns = OrderedDict{String,Vector{Any}}()
    for (name, sym) in zip(names, propertynames(file))
        columns[name] = Any[ismissing(v) ? missing : String(v) for v in getproperty(file, sym)]
    end
    return columns, length(file)
end

"""
    read_saft_csv(path) -> DataTable

`pandas.read_csv(path)`: default NA strings become NaN, and a column whose
values are all numbers (or NaN) becomes numeric.
"""
function read_saft_csv(path)
    columns, n = _csv_strings(path)
    data = Dict{String,Vector{Any}}()
    numeric = Set{String}()
    for (name, values) in columns
        values = Any[(v === missing || v in PANDAS_NA) ? missing : v for v in values]
        if all(v -> v === missing || occursin(PANDAS_NUMBER, v), values)
            data[name] = Any[v === missing ? NaN : py_float(v) for v in values]
            push!(numeric, name)
        else
            data[name] = values
        end
    end
    return DataTable(collect(keys(columns)), data, numeric, collect(0:n-1))
end

"""
    read_text_csv(path, columns) -> DataTable

The user-file reader (`user_species._read_csv`): every value is text (only
empty fields are NaN), the BOM is dropped from the header, and the table has
`columns` first, then any other columns of the file.
"""
function read_text_csv(path, columns)
    raw, n = _csv_strings(path)
    names = unique(vcat(columns, collect(keys(raw))))
    data = Dict{String,Vector{Any}}(c => (haskey(raw, c) ? raw[c] : Any[missing for _ in 1:n]) for c in names)
    return DataTable(names, data, Set{String}(), collect(0:n-1))
end

empty_table(columns) = DataTable(collect(columns), Dict{String,Vector{Any}}(c => Any[] for c in columns),
                                 Set{String}(), Int[])

# ---------------------------------------------------------------- XML
"""
    XMLNode

An element of `supcrt92.xml` with ElementTree's view of it: `text` is the text
before the first child element (`nothing` if there is none).
"""
mutable struct XMLNode
    tag::String
    attributes::OrderedDict{String,String}
    text::Union{Nothing,String}
    children::Vector{XMLNode}
end
XMLNode(tag::AbstractString; text=nothing) = XMLNode(String(tag), OrderedDict{String,String}(), text, XMLNode[])

function _xml_node(e::EzXML.Node)
    attributes = OrderedDict{String,String}(EzXML.nodename(a) => EzXML.nodecontent(a)
                                            for a in EzXML.attributes(e))
    text = nothing
    children = XMLNode[]
    for child in EzXML.eachnode(e)
        if EzXML.iselement(child)
            push!(children, _xml_node(child))
        elseif isempty(children) && (EzXML.istext(child) || EzXML.iscdata(child))
            text = something(text, "") * EzXML.nodecontent(child)
        end
    end
    return XMLNode(EzXML.nodename(e), attributes, text, children)
end

"""Parse an XML file into `XMLNode`s (the root element)."""
read_xml(path) = _xml_node(EzXML.root(EzXML.readxml(path)))

"""First direct child named `tag` (ElementTree `find(tag)`), or nothing."""
findchild(node::XMLNode, tag::AbstractString) = (i = findfirst(c -> c.tag == tag, node.children);
                                                 i === nothing ? nothing : node.children[i])
"""ElementTree `find('a/b')`: a path of direct children."""
function findpath(node::XMLNode, path::AbstractString)
    for tag in split(path, '/')
        (tag == "." || isempty(tag)) && continue
        node = findchild(node, tag)
        node === nothing && return nothing
    end
    return node
end
"""ElementTree `findtext(tag)`: text of the first child `tag` ('' if it has none), nothing without the child."""
function findtext(node::XMLNode, path::AbstractString)
    child = findpath(node, path)
    return child === nothing ? nothing : something(child.text, "")
end
"""ElementTree `find('.//tag')`: first descendant named `tag` in document order, or nothing."""
function finddescendant(node::XMLNode, tag::AbstractString)
    for child in node.children
        child.tag == tag && return child
        found = finddescendant(child, tag)
        found === nothing || return found
    end
    return nothing
end
"""ElementTree `iter(tag)`: the node itself and all descendants named `tag`, in document order."""
function iterxml(node::XMLNode, tag::AbstractString, out=XMLNode[])
    node.tag == tag && push!(out, node)
    for child in node.children
        iterxml(child, tag, out)
    end
    return out
end
getattribute(node::XMLNode, name, default=nothing) = get(node.attributes, name, default)

# ---------------------------------------------------------------- NASA-7
"""One NASA-7 record: temperature limits [K], coefficients a1..a7 of each range, header elements, phase letter."""
struct NASARecord
    T_low::Float64
    T_mid::Float64
    T_high::Float64
    low::Vector{Float64}
    high::Vector{Float64}
    elements::OrderedDict{String,Float64}
    phase::String
end

# ---------------------------------------------------------------- the container
"""
    Databases

The four databases (`nasa`, `pcsaft`, `epcsaft`, `supcrt`), the file `paths`
and the resolved PROJECT_PCSAFT.csv column names `pc_columns`
(m, sigma, epsilon, epsilonAB, kappaAB, association; `nothing` if absent).
"""
mutable struct Databases
    nasa::OrderedDict{String,NASARecord}
    pcsaft::DataTable
    epcsaft::DataTable
    supcrt::XMLNode
    paths::OrderedDict{String,String}
    pc_columns::OrderedDict{String,Union{Nothing,String}}
end

"""
    load_databases(data_dir; paths=nothing) -> Databases

Read all databases from `data_dir` (default file names) or explicit `paths`
(a dict with the keys of `DATABASE_FILES`). Python: `Databases.load`.
"""
function load_databases(data_dir=nothing; paths=nothing)
    if paths === nothing
        paths = OrderedDict(key => joinpath(String(data_dir), name) for (key, name) in DATABASE_FILES)
    end
    paths = OrderedDict{String,String}(key => String(paths[key]) for key in keys(DATABASE_FILES))
    for (key, path) in paths
        isfile(path) || throw(ArgumentError("$key database not found: $path"))
    end
    pcsaft = read_saft_csv(paths["pcsaft"])
    nasa = read_nasa_database(paths["nasa"])
    epcsaft = read_saft_csv(paths["epcsaft"])
    supcrt = read_xml(paths["supcrt"])
    pc_columns = OrderedDict{String,Union{Nothing,String}}(
        "m" => "m", "sigma" => "sigma_A", "epsilon" => "u_over_kB_K",
        "epsilonAB" => first_existing_column(pcsaft, ["epsilonAB_over_kB_K", "epsilonAB_K",
                                                      "epsilonAB/kB", "epsilonAB"]),
        "kappaAB" => first_existing_column(pcsaft, ["kappaAB", "association_volume_kappaAB"]),
        "association" => first_existing_column(pcsaft, ["Association_scheme", "association_scheme",
                                                        "Association"]),
    )
    return Databases(nasa, pcsaft, epcsaft, supcrt, paths, pc_columns)
end

Base.copy(db::Databases) = deepcopy(db)

const _SUBSCRIPTS = ('₂' => '2', '₃' => '3', '₄' => '4', '₅' => '5', '₆' => '6', '₇' => '7', '₈' => '8', '₉' => '9')

"""
    normalize_name(value)

Case-folded name with subscript digits replaced, for matching CSV entries
('' for NaN/missing).
"""
function normalize_name(value)
    isna(value) && return ""
    return replace(casefold(strip(py_str(value))), _SUBSCRIPTS...)
end

function first_existing_column(table::DataTable, candidates)
    normalized = Dict{String,String}()
    for col in table.columns
        normalized[normalize_name(col)] = col
    end
    for candidate in candidates
        key = normalize_name(candidate)
        haskey(normalized, key) && return normalized[key]
    end
    return nothing
end

# ---------------------------------------------------------------- NASA-7 parser
function nasa_float(text)
    text = replace(strip(py_str(text)), 'D' => 'E', 'd' => 'e')
    isempty(text) && return NaN
    return py_float(text)
end

"""
    parse_nasa_header_elements(header) -> OrderedDict{String,Float64}

Elements in the four 5-column fields of a CHEMKIN header (columns 25-44).
"""
function parse_nasa_header_elements(header::AbstractString)
    elements = OrderedDict{String,Float64}()
    for offset in 24:5:39
        field = py_slice(header, offset, offset + 5)
        symbol = strip(py_slice(field, 0, 2))
        symbol = isempty(symbol) ? "" : uppercase(first(symbol)) * lowercase(symbol[nextind(symbol, 1):end])
        isempty(symbol) && continue
        count = py_tryfloat(String(strip(py_slice(field, 2, 5))))
        count === nothing && continue
        if count != 0
            elements[symbol] = get(elements, symbol, 0.0) + count
        end
    end
    return elements
end

"""Parse 4-line NASA-7 records into an OrderedDict name => NASARecord; error if there are none."""
function read_nasa_database(path)
    nasa_db = read_nasa_lines(py_splitlines(read_text(path); keepends=true))
    isempty(nasa_db) && throw(ArgumentError("No NASA-7 entries found"))
    return nasa_db
end

"""
    read_nasa_lines(lines) -> OrderedDict{String,NASARecord}

The NASA-7 records found in `lines` (strict CHEMKIN columns, markers 1-4 in
column 80); empty if there are none. A later record with the same name
replaces the earlier one in its position (as a Python dict).
"""
function read_nasa_lines(lines::AbstractVector{<:AbstractString})
    nasa_db = OrderedDict{String,NASARecord}()
    marker(line, k) = length(line) >= 80 && py_slice(line, 79, 80) == string(k)
    i = 1
    while i <= length(lines)
        line1 = rstrip(lines[i], '\n')
        if !marker(line1, 1)
            i += 1
            continue
        end
        i + 3 > length(lines) && break
        line2, line3, line4 = (rstrip(lines[i+k], '\n') for k in 1:3)
        if !(marker(line2, 2) && marker(line3, 3) && marker(line4, 4))
            i += 1
            continue
        end
        species_name = String(strip(py_slice(line1, 0, 18)))
        T_low, T_high, T_mid = try
            nasa_float(py_slice(line1, 45, 55)), nasa_float(py_slice(line1, 55, 65)), nasa_float(py_slice(line1, 65, 73))
        catch
            throw(ArgumentError("Could not parse NASA temperature limits for species '$species_name'." *
                                "\nHeader:\n$line1"))
        end
        high, low = try
            c2 = [nasa_float(py_slice(line2, j, j + 15)) for j in 0:15:60]
            c3 = [nasa_float(py_slice(line3, j, j + 15)) for j in 0:15:60]
            c4 = [nasa_float(py_slice(line4, j, j + 15)) for j in 0:15:45]
            [c2[1], c2[2], c2[3], c2[4], c2[5], c3[1], c3[2]], [c3[3], c3[4], c3[5], c4[1], c4[2], c4[3], c4[4]]
        catch
            throw(ArgumentError("Could not parse NASA coefficients for species '$species_name'."))
        end
        nasa_db[species_name] = NASARecord(T_low, T_mid, T_high, low, high, parse_nasa_header_elements(line1),
                                           uppercase(strip(py_slice(line1, 44, 45))))
        i += 4
    end
    return nasa_db
end

function nasa7_coefficients(nasa_db, species_name, T)
    haskey(nasa_db, species_name) || throw(KeyError("NASA species not found: $species_name"))
    data = nasa_db[species_name]
    if T < data.T_low || T > data.T_high
        throw(ArgumentError(@sprintf("Temperature %.2f K is outside the NASA range for %s: %.2f–%.2f K.",
                                     T, species_name, data.T_low, data.T_high)))
    end
    T <= data.T_mid && return data.low, "low"
    return data.high, "high"
end

"""
    nasa7_properties(nasa_db, species_name, T) -> (cp [J/mol/K], h [J/mol], s [J/mol/K], g = h - T s [J/mol])

From the NASA-7 polynomials at T [K].
"""
function nasa7_properties(nasa_db, species_name, T)
    T = Float64(T)
    T <= 0.0 && throw(ArgumentError("NASA temperature must be positive."))
    a, _ = nasa7_coefficients(nasa_db, species_name, T)
    a1, a2, a3, a4, a5, a6, a7 = a
    T2, T3, T4 = T^2, T^3, T^4
    cp_R = a1 + a2 * T + a3 * T2 + a4 * T3 + a5 * T4
    h_RT = a1 + a2 * T / 2.0 + a3 * T2 / 3.0 + a4 * T3 / 4.0 + a5 * T4 / 5.0 + a6 / T
    s_R = a1 * log(T) + a2 * T + a3 * T2 / 2.0 + a4 * T3 / 3.0 + a5 * T4 / 4.0 + a7
    cp = R * cp_R
    h = R * T * h_RT
    entropy = R * s_R
    g = h - T * entropy
    return cp, h, entropy, g
end
