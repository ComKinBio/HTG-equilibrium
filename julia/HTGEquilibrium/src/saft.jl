# Selection of PC-SAFT / ePC-SAFT parameter sets and SUPCRT HKF records
# (Python: htg/saft.py). Only the parameter data and its selection; building
# the Clapeyron models is part of phases.jl (a later step).
#
# Molecules: the first complete ePC-SAFT set (PROJECT_EPCSAFT.csv) is used,
# otherwise the first complete PC-SAFT set (PROJECT_PCSAFT.csv, model priority
# 'PC-SAFT' then 'GC-PPC-SAFT'). Water's ePC-SAFT sigma is temperature dependent:
# sigma(T) = 2.7927 + 10.11 exp(-0.01775 T) - 1.417 exp(-0.01146 T) Angstrom, T in K.
#
# Ions: the first complete charged ePC-SAFT set, plus the first complete aqueous
# HKF record in supcrt92.xml.
#
# Parameter sets are returned as OrderedDict{String,Any} with Python's keys
# ('Species', 'Formula', 'Source', 'Database row', 'Model', 'RefID', 'm', ...).

const EPC_PARAMETER_ALIASES = OrderedDict(
    "m" => ["mseg", "m", "segment"],
    "sigma" => ["sigma", "sigma_a", "sigma (a)", "sigma (å)", "σ"],
    "epsilon" => ["u/kb", "epsilon/k", "epsilon/kb", "epsilon_k", "u_over_kb_k"],
    "epsilonAB" => ["epsilonab/kb", "association_energy_epsilonab/k", "epsilonab_over_kb_k", "epsilonab_k"],
    "kappaAB" => ["kappaab", "association_volume_kappaab"],
    "association_scheme" => ["association_scheme"],
)
const ION_EPC_PARAMETER_ALIASES = OrderedDict(
    "m" => ["mseg", "m", "segment"],
    "sigma" => ["sigma", "sigma_a", "sigma (a)", "sigma (å)"],
    "epsilon" => ["u/kb", "epsilon/k", "epsilon/kb", "epsilon_k"],
    "dielectric_constant" => ["dielectric_constant"],
)
const PCSAFT_MODEL_PRIORITY = ["pc-saft", "gc-ppc-saft"]
const HKF_REQUIRED_FIELDS = ["Gf", "Hf", "Sr", "a1", "a2", "a3", "a4", "c1", "c2", "wref"]

const ParameterSet = OrderedDict{String,Any}

"""File a database row came from: rows of the user layer carry a 'source_file' column."""
_source(t::DataTable, i, default) = (v = cell(t, "source_file", i); isna(v) ? default : py_str(v))

"""Temperature-dependent ePC-SAFT segment diameter of water [Angstrom] at T [K]."""
water_sigma(T) = Float64(2.7927 + 10.11 * exp(-0.01775 * T) - 1.417 * exp(-0.01146 * T))

_sort_by_index(t::DataTable, rows) = sort(rows; by=i -> t.index[i])
_text_eq(v, s) = !isna(v) && casefold(strip(py_str(v))) == s   # astype(str).str.strip().str.casefold().eq(s)

function _parameter_rows(t::DataTable, rows)
    hascolumn(t, "record_type") || return rows
    selected = [i for i in rows if _text_eq(cell(t, "record_type", i), "species_parameter")]
    return isempty(selected) ? rows : selected
end

# ---------------------------------------------------------------- molecules
function _evaluate_epc_value(species_name, quantity, raw, T)
    isna(raw) && return NaN
    quantity == "association_scheme" && return String(strip(py_str(raw)))
    v = py_tryfloat(raw)
    v === nothing || return v
    text = strip(py_str(raw))
    if species_name == "water" && quantity == "sigma" && startswith(casefold(text), "t-dependent")
        return water_sigma(T)
    end
    throw(ArgumentError("Cannot evaluate ePC-SAFT parameter '$quantity' for '$species_name'.\n" *
                        "Raw value: $(py_repr(raw))"))
end

function _extract_epc_quantity(t::DataTable, rows, species_name, quantity, T; required=true)
    aliases = Set(normalize_name(a) for a in EPC_PARAMETER_ALIASES[quantity])
    matches = [i for i in rows if normalize_name(cell(t, "Parameter", i)) in aliases]
    isempty(matches) && return required ? nothing : NaN
    try
        return _evaluate_epc_value(species_name, quantity, cell(t, "Value", matches[1]), T)
    catch
        return required ? nothing : NaN
    end
end

"""Normalized names that identify a species in the SAFT CSV Species column."""
function species_aliases(species_name, aliases, formula)
    names = Set(normalize_name(v) for v in (aliases === nothing ? [species_name] : aliases))
    formula === nothing || push!(names, normalize_name(formula))
    return names
end

"""Groups of `rows` with equal values in `columns` (NaN equal to NaN), in order of first appearance."""
function _groups(t::DataTable, rows, columns)
    isempty(columns) && return [rows]
    keys_ = Any[]
    groups = Vector{Int}[]
    for i in rows
        key = Tuple(isna(cell(t, c, i)) ? missing : cell(t, c, i) for c in columns)
        k = findfirst(x -> isequal(x, key), keys_)
        if k === nothing
            push!(keys_, key)
            push!(groups, [i])
        else
            push!(groups[k], i)
        end
    end
    return groups
end

"""
    first_complete_epcsaft_set(epcsaft, species_name, formula, aliases, T)

First (Model, RefID) group of ePC-SAFT rows of the species with finite m,
sigma, epsilon at T [K]; `nothing` if there is none.
"""
function first_complete_epcsaft_set(epcsaft::DataTable, species_name, formula, aliases, T)
    names = species_aliases(species_name, aliases, formula)
    rows = [i for i in 1:nrows(epcsaft) if normalize_name(cell(epcsaft, "Species", i)) in names]
    isempty(rows) && return nothing
    rows = _sort_by_index(epcsaft, _parameter_rows(epcsaft, rows))
    grouping = [c for c in ("Model", "RefID") if hascolumn(epcsaft, c)]
    for group in _groups(epcsaft, rows, grouping)
        m_i = _extract_epc_quantity(epcsaft, group, species_name, "m", T)
        sigma_i = _extract_epc_quantity(epcsaft, group, species_name, "sigma", T)
        epsilon_i = _extract_epc_quantity(epcsaft, group, species_name, "epsilon", T)
        (m_i === nothing || sigma_i === nothing || epsilon_i === nothing) && continue
        all(isfinite, (m_i, sigma_i, epsilon_i)) || continue
        epsilonAB_i = _extract_epc_quantity(epcsaft, group, species_name, "epsilonAB", T; required=false)
        kappaAB_i = _extract_epc_quantity(epcsaft, group, species_name, "kappaAB", T; required=false)
        association_i = _extract_epc_quantity(epcsaft, group, species_name, "association_scheme", T;
                                              required=false)
        formula === nothing && throw(KeyError(species_name))
        first_row = group[1]
        refid = cell(epcsaft, "RefID", first_row)
        return ParameterSet(
            "Species" => species_name, "Formula" => formula,
            "Source" => _source(epcsaft, first_row, "PROJECT_EPCSAFT.csv"),
            "Database row" => minimum(epcsaft.index[i] for i in group),
            "Model" => hascolumn(epcsaft, "Model") ? py_str(cell(epcsaft, "Model", first_row)) : "ePC-SAFT",
            "RefID" => (hascolumn(epcsaft, "RefID") && !isna(refid)) ? py_str(refid) : "",
            "m" => Float64(m_i), "sigma" => Float64(sigma_i), "epsilon" => Float64(epsilon_i),
            "epsilonAB" => isna(epsilonAB_i) ? NaN : Float64(epsilonAB_i),
            "kappaAB" => isna(kappaAB_i) ? NaN : Float64(kappaAB_i),
            "association_scheme" => isna(association_i) ? "-" : String(strip(py_str(association_i))),
        )
    end
    return nothing
end

function _optional_numeric(t::DataTable, i, column)
    column === nothing && return NaN
    value = cell(t, column, i)
    isna(value) && return NaN
    v = py_tryfloat(value)
    return v === nothing ? NaN : v
end

"""
    first_complete_pcsaft_set(databases, species_name, formula)

First positive, finite PC-SAFT row matching `formula` (model priority PC-SAFT,
then GC-PPC-SAFT); `nothing` if there is none.
"""
function first_complete_pcsaft_set(databases::Databases, species_name, formula)
    formula === nothing && throw(KeyError("No PC-SAFT formula mapping for '$species_name'."))
    pcsaft, cols = databases.pcsaft, databases.pc_columns
    key = normalize_name(formula)
    rows = [i for i in 1:nrows(pcsaft) if normalize_name(cell(pcsaft, "Formula", i)) == key]
    isempty(rows) && return nothing
    rows = _sort_by_index(pcsaft, rows)
    value(i, k) = (c = cols[k]; (c === nothing || !hascolumn(pcsaft, c)) ? nothing :
                   (v = cell(pcsaft, c, i); v === missing ? NaN : py_tryfloat(v)))
    for preferred in PCSAFT_MODEL_PRIORITY
        for i in rows
            _text_eq(cell(pcsaft, "Model", i), preferred) || continue
            m_i, sigma_i, epsilon_i = value(i, "m"), value(i, "sigma"), value(i, "epsilon")
            (m_i === nothing || sigma_i === nothing || epsilon_i === nothing) && continue
            all(isfinite, (m_i, sigma_i, epsilon_i)) || continue
            (m_i <= 0.0 || sigma_i <= 0.0 || epsilon_i <= 0.0) && continue
            association_i = "-"
            assoc = cols["association"]
            if assoc !== nothing && !isna(cell(pcsaft, assoc, i))
                association_i = String(strip(py_str(cell(pcsaft, assoc, i))))
            end
            refid = cell(pcsaft, "RefID", i)
            return ParameterSet(
                "Species" => species_name, "Formula" => formula,
                "Source" => _source(pcsaft, i, "PROJECT_PCSAFT.csv"),
                "Database row" => pcsaft.index[i], "Model" => py_str(cell(pcsaft, "Model", i)),
                "RefID" => (hascolumn(pcsaft, "RefID") && !isna(refid)) ? py_str(refid) : "",
                "m" => m_i, "sigma" => sigma_i, "epsilon" => epsilon_i,
                "epsilonAB" => _optional_numeric(pcsaft, i, cols["epsilonAB"]),
                "kappaAB" => _optional_numeric(pcsaft, i, cols["kappaAB"]),
                "association_scheme" => association_i,
            )
        end
    end
    return nothing
end

# ---------------------------------------------------------------- ions
"""Ion name for matching: case-folded, spaces removed, ₂ ₃ ⁺ ⁻ replaced ('' for NaN)."""
function normalize_ion_name(value)
    isna(value) && return ""
    return replace(casefold(strip(py_str(value))), ' ' => "", '₂' => '2', '₃' => '3', '⁺' => '+', '⁻' => '-')
end

function _numeric_epc_value(value)
    isna(value) && return NaN
    v = py_tryfloat(value)
    v === nothing && throw(ArgumentError("Could not convert ePC-SAFT value $(py_repr(value)) to a number."))
    return v
end

function _first_epc_parameter(t::DataTable, rows, quantity; required=true)
    for alias in (casefold(strip(x)) for x in ION_EPC_PARAMETER_ALIASES[quantity])
        for i in rows
            _text_eq(cell(t, "Parameter", i), alias) && return _numeric_epc_value(cell(t, "Value", i))
        end
    end
    if required
        available = unique([py_str(cell(t, "Parameter", i)) for i in rows if !isna(cell(t, "Parameter", i))])
        throw(ArgumentError("Required ePC-SAFT parameter '$quantity' not found.\n" *
                            "Available parameters: [$(join(py_repr.(available), ", "))]"))
    end
    return NaN
end

"""
    first_complete_ion_epcsaft_set(epcsaft, ion_name, saft_formula, charge)

First (RefID, Model) set of charged ePC-SAFT rows of the ion with positive,
finite m, sigma, epsilon; `nothing` if there is none.
"""
function first_complete_ion_epcsaft_set(epcsaft::DataTable, ion_name, saft_formula, charge)
    aliases = Set([normalize_ion_name(ion_name), normalize_ion_name(saft_formula)])
    rows = [i for i in 1:nrows(epcsaft) if normalize_ion_name(cell(epcsaft, "Species", i)) in aliases]
    isempty(rows) && return nothing
    rows = _parameter_rows(epcsaft, _sort_by_index(epcsaft, rows))
    grouping = [c for c in ("RefID", "Model") if hascolumn(epcsaft, c)]
    for parameter_set in _groups(epcsaft, rows, grouping)
        m_i, sigma_i, epsilon_i = try
            (_first_epc_parameter(epcsaft, parameter_set, "m"), _first_epc_parameter(epcsaft, parameter_set, "sigma"),
             _first_epc_parameter(epcsaft, parameter_set, "epsilon"))
        catch err
            err isa ArgumentError || rethrow()
            continue
        end
        all(isfinite, (m_i, sigma_i, epsilon_i)) || continue
        (m_i <= 0.0 || sigma_i <= 0.0 || epsilon_i <= 0.0) && continue
        dielectric_i = _first_epc_parameter(epcsaft, parameter_set, "dielectric_constant"; required=false)
        first_row = parameter_set[1]
        charge_database = NaN
        if hascolumn(epcsaft, "Charge")
            values = [cell(epcsaft, "Charge", i) for i in parameter_set if !isna(cell(epcsaft, "Charge", i))]
            isempty(values) || (charge_database = py_float(values[1]))
        end
        model, refid = cell(epcsaft, "Model", first_row), cell(epcsaft, "RefID", first_row)
        return ParameterSet(
            "Species" => ion_name, "Source" => _source(epcsaft, first_row, "PROJECT_EPCSAFT.csv"),
            "Database row" => minimum(epcsaft.index[i] for i in parameter_set),
            "Model" => isna(model) ? "" : py_str(model), "RefID" => isna(refid) ? "" : py_str(refid),
            "charge" => charge, "charge_database" => charge_database,
            "m" => m_i, "sigma" => sigma_i, "epsilon" => epsilon_i, "dielectric_constant" => dielectric_i,
        )
    end
    return nothing
end

function supcrt_species_name(node::XMLNode)
    name = findchild(node, "Name")
    (name === nothing || name.text === nothing) && return ""
    return String(strip(name.text))
end

"""Charge and the aqueous HKF fields of one SUPCRT record; `nothing` if incomplete."""
function read_hkf_candidate(node::XMLNode)
    charge_element = findchild(node, "Charge")
    (charge_element === nothing || charge_element.text === nothing) && return nothing
    charge = py_tryfloat(strip(charge_element.text))
    charge === nothing && return nothing
    hkf = findpath(node, "Thermo/HKF")
    hkf === nothing && return nothing
    result = OrderedDict{String,Any}("Charge" => charge)
    for name in HKF_REQUIRED_FIELDS
        child = findchild(hkf, name)
        (child === nothing || child.text === nothing || isempty(strip(child.text))) && return nothing
        value = py_tryfloat(strip(child.text))
        value === nothing && return nothing
        result[name] = value
    end
    return result
end

"""
    first_complete_hkf_set(supcrt_root, ion_name, supcrt_aliases)

First SUPCRT record named like one of `supcrt_aliases` with a complete aqueous
HKF record: `Species`, `SUPCRT name`, `Candidate` (0-based), `Charge`, the HKF
fields, `Formula`, `MolarMass`, `Type`, `Date`, `References`; `nothing` if none.
"""
function first_complete_hkf_set(supcrt_root::XMLNode, ion_name, supcrt_aliases)
    aliases = Set(normalize_ion_name(a) for a in supcrt_aliases)
    candidates = [node for node in iterxml(supcrt_root, "Species")
                  if normalize_ion_name(supcrt_species_name(node)) in aliases]
    for (number, candidate) in enumerate(candidates)
        values = read_hkf_candidate(candidate)
        values === nothing && continue
        result = ParameterSet("Species" => ion_name, "SUPCRT name" => supcrt_species_name(candidate),
                              "Candidate" => number - 1)
        merge!(result, values)
        for tag in ("Formula", "MolarMass", "Type", "Date", "References")
            child = finddescendant(candidate, tag)
            result[tag] = (child !== nothing && child.text !== nothing) ? String(strip(child.text)) : ""
        end
        return result
    end
    return nothing
end
