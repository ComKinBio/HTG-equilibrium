# Species registry: which species exist, their composition and database
# identifiers (Python: htg/species.py).
#
# The registry starts from data/species.toml and is completed by scanning the
# databases (`discover!`). It never invents data: a species is "ready" only
# when every required record exists. `registry.catalogue` lists every
# candidate, ready or not, with the reason.
#
# Registry entries are OrderedDict{String,Any} with Python's keys (formula,
# nasa, elements; solids from supcrt92.xml also adapter, supcrt_mineral,
# source, references, date, hydrate; ions formula, charge, supcrt, elements).
# As in Python, a molecule discovered as ready is the same entry object in
# `molecular` and `organic`. The registry is filled while the model loads and
# read only afterwards.

const ELEMENT_PATTERN = r"([A-Z][a-z]?)\(([-+\d.eE]+)\)"

const RegistryEntry = OrderedDict{String,Any}

"""
    formula_elements(formula) -> OrderedDict or nothing

Elemental composition of a plain formula like 'C2H6O'; `nothing` if it cannot
be parsed (Python: species.formula_elements).
"""
function formula_elements(formula)
    text = strip(py_str(formula))
    occursin(r"^(?:[A-Z][a-z]?\d*)+\z", text) || return nothing
    result = Elements()
    for m in eachmatch(r"([A-Z][a-z]?)(\d*)", text)
        result[m[1]] = get(result, m[1], 0.0) + (isempty(m[2]) ? 1 : parse(Int, m[2]))
    end
    return result
end

"""Molar mass [g/mol] of a plain formula (Python: species.molecular_weight)."""
function molecular_weight(formula, atomic_weights=ATOMIC_WEIGHTS)
    tokens = collect(eachmatch(r"([A-Z][a-z]?)(\d*)", strip(py_str(formula))))
    isempty(tokens) && throw(ArgumentError("Could not parse molecular formula: $formula"))
    total = 0.0
    for m in tokens
        element, number = m[1], m[2]
        haskey(atomic_weights, element) ||
            throw(KeyError("Atomic weight missing for '$element' in formula '$formula'."))
        total += atomic_weights[element] * (isempty(number) ? 1 : parse(Int, number))
    end
    return Float64(total)
end

function _ion_identity(name, charge)
    text = replace(replace(py_str(name), r"\s+" => ""), "^" => "")
    # The explicit charge avoids treating the 2 in VO2+ as a divalent charge.
    sign = charge > 0 ? "+" : "-"
    n = abs(Int(trunc(charge)))
    suffixes = ["$n$sign", "$sign$n", sign^n]
    if uppercase(text) in ("NA+", "K+", "LI+", "RB+", "CS+", "CL-", "F-", "BR-", "I-")
        return casefold(rstrip(text, ('+', '-'))), Int(trunc(charge))
    end
    for suffix in suffixes
        if endswith(text, suffix)
            return casefold(first(text, length(text) - length(suffix))), Int(trunc(charge))
        end
    end
    return casefold(rstrip(text, ('+', '-'))), Int(trunc(charge))
end

const NO_SAFT_REASON = "NASA available; no resolved SAFT mapping. Add SAFT parameters, or a " *
                       "[nasa_names] entry in species.toml (see docs/adding_species.md)."

"""One row of the species catalogue (Python: a row of `Registry.catalogue`). `nothing` is Python's None."""
struct CatalogueRow
    Species::Union{Nothing,String}
    Formula::Union{Nothing,String}
    Phase::String
    Ready::Bool
    Reason::String
    NASA::Union{Nothing,String}
end
CatalogueRow(; Species, Formula, Phase, Ready, Reason, NASA) =
    CatalogueRow(Species, Formula, Phase, Ready, Reason, NASA)

"""
    Registry(config)

All species known to the model, from a parsed species.toml (`config`, an
OrderedDict as from `parsefile_ordered`):

    molecular, organic  fluid molecules {name => {formula, nasa, elements}}
    solid               pure solids {name => {formula, nasa | supcrt_mineral, elements, ...}}
    ion                 aqueous ions {name => {formula, charge, supcrt, elements}}

plus the identifiers used to look species up in the SAFT and SUPCRT databases,
the aqueous reference components, binary interactions, the `catalogue` and the
`ion_epcsaft_table` / `ion_hkf_table` rows (filled by `discover!`).
"""
mutable struct Registry
    config::OrderedDict{String,Any}
    molecular::OrderedDict{String,RegistryEntry}
    organic::OrderedDict{String,RegistryEntry}
    solid::OrderedDict{String,RegistryEntry}
    ion::OrderedDict{String,RegistryEntry}
    saft_formula::OrderedDict{String,String}
    saft_aliases::OrderedDict{String,Vector{String}}
    ion_charge::OrderedDict{String,Int}
    ion_saft_formula::OrderedDict{String,String}
    supcrt_ion_aliases::OrderedDict{String,Vector{Union{Nothing,String}}}
    aqueous_reference_neutrals::Vector{String}
    aqueous_reference_ions::Vector{String}
    probe_weights::OrderedDict{String,Float64}
    binary_interactions::Vector{Tuple{Tuple{String,String},Float64}}
    catalogue::Vector{CatalogueRow}
    ion_epcsaft_table::Vector{ParameterSet}
    ion_hkf_table::Vector{ParameterSet}
end

_elements_of(d) = Elements(String(k) => Float64(v) for (k, v) in d)

function Registry(config::AbstractDict)
    config = OrderedDict{String,Any}(config)
    reg = Registry(config, OrderedDict{String,RegistryEntry}(), OrderedDict{String,RegistryEntry}(),
                   OrderedDict{String,RegistryEntry}(), OrderedDict{String,RegistryEntry}(),
                   OrderedDict{String,String}(), OrderedDict{String,Vector{String}}(),
                   OrderedDict{String,Int}(), OrderedDict{String,String}(),
                   OrderedDict{String,Vector{Union{Nothing,String}}}(), String[], String[],
                   OrderedDict{String,Float64}(), Tuple{Tuple{String,String},Float64}[],
                   CatalogueRow[], ParameterSet[], ParameterSet[])
    for group in ("molecular", "organic")
        table = getfield(reg, Symbol(group))
        for (name, entry) in get(config, group, Dict())
            table[name] = RegistryEntry("formula" => entry["formula"], "nasa" => entry["nasa"],
                                        "elements" => _elements_of(entry["elements"]))
            reg.saft_formula[name] = entry["formula"]
            reg.saft_aliases[name] = String.(entry["saft_aliases"])
        end
    end
    for (name, entry) in get(config, "solid", Dict())
        reg.solid[name] = RegistryEntry("formula" => entry["formula"], "nasa" => entry["nasa"],
                                        "elements" => _elements_of(entry["elements"]))
    end
    for (name, entry) in get(config, "ion", Dict())
        reg.ion[name] = RegistryEntry("formula" => entry["formula"], "charge" => entry["charge"],
                                      "supcrt" => entry["supcrt"], "elements" => _elements_of(entry["elements"]))
        reg.ion_charge[name] = entry["charge"]
        reg.ion_saft_formula[name] = entry["saft_formula"]
        reg.supcrt_ion_aliases[name] = collect(Union{Nothing,String}, entry["supcrt_aliases"])
    end
    reference = config["aqueous_reference"]
    append!(reg.aqueous_reference_neutrals, reference["neutrals"])
    append!(reg.aqueous_reference_ions, reference["ions"])
    for (k, v) in reference["probe_weights"]
        reg.probe_weights[k] = Float64(v)
    end
    for item in get(config, "binary_interaction", [])
        a, b = item["components"]
        push!(reg.binary_interactions, ((String(a), String(b)), Float64(item["kij"])))
    end
    return reg
end

"""Registry from a species.toml file; `nasa_names` replaces its [nasa_names] table."""
function registry_from_file(path; nasa_names=nothing)
    config = parsefile_ordered(path)
    nasa_names === nothing || (config["nasa_names"] = OrderedDict{String,Any}(nasa_names))
    return Registry(config)
end

Base.copy(reg::Registry) = deepcopy(reg)

# ---------------------------------------------------------------- lookups
"""Fluid molecules: `molecular` and `organic` (organic entries win), as Python's {**molecular, **organic}."""
all_molecular(reg::Registry) = merge(reg.molecular, reg.organic)

"""NASA record name of a molecule or solid."""
function nasa_name(reg::Registry, name)
    for table in (reg.solid, reg.organic, reg.molecular)   # later Python dicts win
        haskey(table, name) && return table[name]["nasa"]
    end
    throw(KeyError("No NASA mapping for $name"))
end

"""Registry entry of `name` in `phase`; error with the available names otherwise."""
function entry(reg::Registry, phase, name)
    registry = phase == "solid" ? reg.solid :
               (phase == "aqueous" && haskey(reg.ion, name)) ? reg.ion : all_molecular(reg)
    haskey(registry, name) ||
        throw(ArgumentError("$phase: '$name' is not registered. Available species: $(join(keys(registry), ", "))"))
    return registry[name]
end

any_entry(reg::Registry, name) = merge(all_molecular(reg), reg.solid, reg.ion)[name]

"""Catalogue rows: ready species only, or every candidate with `include_unavailable`."""
available(reg::Registry; include_unavailable=false) =
    include_unavailable ? copy(reg.catalogue) : filter(r -> r.Ready, reg.catalogue)

# ---------------------------------------------------------------- SAFT parameters
"""First complete ePC-SAFT set, else PC-SAFT set, for molecule `name` at T [K]; `nothing` if none."""
function saft_parameters(reg::Registry, databases::Databases, name, T)
    formula = get(reg.saft_formula, name, nothing)
    data = first_complete_epcsaft_set(databases.epcsaft, name, formula, get(reg.saft_aliases, name, nothing), T)
    return data === nothing ? first_complete_pcsaft_set(databases, name, formula) : data
end

"""
    resolve_parameters(registry, databases, T, species; atomic_weights=ATOMIC_WEIGHTS)

SAFT parameters of every molecule in `species` at T [K] (OrderedDict name =>
parameter set, with `Mw` added); error when a set is incomplete. Water without
an association scheme gets '2B'.
"""
function resolve_parameters(reg::Registry, databases::Databases, T, species; atomic_weights=ATOMIC_WEIGHTS)
    resolved = OrderedDict{String,ParameterSet}()
    molecules = all_molecular(reg)
    for name in unique(species)
        haskey(molecules, name) || throw(ArgumentError("$name: not found in the molecular database registry"))
        data = saft_parameters(reg, databases, name, T)
        data === nothing && throw(ArgumentError("$name: NASA data alone is insufficient for PC-SAFT; " *
                                                "complete SAFT parameters/mapping required"))
        data["Mw"] = sum(atomic_weights[e] * v for (e, v) in molecules[name]["elements"])
        if name == "water" && get(data, "association_scheme", nothing) in ("-", "", nothing)
            data["association_scheme"] = "2B"
        end
        resolved[name] = data
    end
    return resolved
end

# ---------------------------------------------------------------- discovery
_truthy(x) = !(x === nothing || x === missing || x == "")
_has(db, key) = key isa AbstractString && haskey(db, key)
_charge_value(v) = isna(v) ? 0.0 : v   # pandas fillna(0)

"""
    discover!(registry, databases, T)

Add every species with complete database records and build the catalogue.
T [K] is used only to evaluate temperature-dependent SAFT parameters while
checking completeness. Python: `Registry.discover`.
"""
function discover!(reg::Registry, databases::Databases, T)
    AW = ATOMIC_WEIGHTS
    nasa_db = databases.nasa
    overrides = Dict(casefold(strip(py_str(k))) => v for (k, v) in get(reg.config, "nasa_names", Dict()))
    original_molecules = all_molecular(reg)
    nasa_by_formula = OrderedDict{Vector{Pair{String,Float64}},Vector{String}}()
    for (name, data) in nasa_db
        key = sort(collect(data.elements); by=first)
        isempty(key) || push!(get!(nasa_by_formula, key, String[]), name)
    end
    pcsaft, epcsaft = databases.pcsaft, databases.epcsaft
    pc = rows_where(pcsaft, "record_type", "pure_component")
    epc = [i for i in rows_where(epcsaft, "record_type", "species_parameter")
           if isequal(_charge_value(cell(epcsaft, "Charge", i)), 0.0)]
    records = Tuple{String,String}[]
    seen = Set{String}()
    for (frame, rows) in ((pcsaft, pc), (epcsaft, epc))
        for i in rows
            (isna(cell(frame, "Species", i)) || isna(cell(frame, "Formula", i))) && continue
            name = casefold(strip(py_str(cell(frame, "Species", i))))
            if !(name in seen)
                push!(seen, name)
                push!(records, (name, String(strip(py_str(cell(frame, "Formula", i))))))
            end
        end
    end

    # Fluid molecules: SAFT record + unambiguous NASA record.
    catalogue = CatalogueRow[]
    for (name, formula) in records
        elements = formula_elements(formula)
        override = get(overrides, name, nothing)
        mapped = _truthy(override) ? override :
                 get(get(original_molecules, name, RegistryEntry()), "nasa", nothing)
        if !_truthy(mapped)
            exact = nothing
            for key in keys(nasa_db)
                if casefold(key) in (name, casefold(formula))
                    exact = key
                    break
                end
            end
            if _truthy(exact)
                mapped = exact
            elseif elements !== nothing && !isempty(elements)
                candidates = get(nasa_by_formula, sort(collect(elements); by=first), String[])
                same_formula = [other for (other, f) in records
                                if (fe = formula_elements(f); fe !== nothing && Dict(fe) == Dict(elements))]
                if length(candidates) == 1 && length(same_formula) == 1
                    mapped = candidates[1]
                end
            end
        end
        if _truthy(override) && !_has(nasa_db, mapped)
            throw(ArgumentError("NASA override for $name not found: $(py_str(mapped))"))
        end
        if _has(nasa_db, mapped) && !isempty(nasa_db[mapped].elements)
            elements = nasa_db[mapped].elements
        end
        mapped_name = mapped === nothing ? nothing : String(mapped)
        if elements === nothing || isempty(elements)
            push!(catalogue, CatalogueRow(name, formula, "molecular", false, "Elemental composition unresolved",
                                          mapped_name))
            continue
        end
        reg.saft_formula[name] = formula
        reg.saft_aliases[name] = [name, formula]
        data = saft_parameters(reg, databases, name, T)
        ready = _has(nasa_db, mapped) && data !== nothing && all(e -> haskey(AW, e), keys(elements))
        reason = ready ? "ready" : !_has(nasa_db, mapped) ? "NASA mapping missing/ambiguous" :
                 data === nothing ? "Complete SAFT parameters missing" : "Atomic weight missing"
        push!(catalogue, CatalogueRow(name, formula, "molecular", ready, reason, mapped_name))
        if ready
            e = RegistryEntry("formula" => formula, "nasa" => mapped_name, "elements" => elements)
            reg.molecular[name] = e
            reg.organic[name] = e
        end
    end

    # Every NASA record stays usable as a dry-feed species; solid records are pure solids.
    covered = Set(get(e, "nasa", nothing) for e in values(reg.molecular))
    for (name, data) in nasa_db
        elements = data.elements
        (isempty(elements) || !all(e -> haskey(AW, e), keys(elements))) && continue
        key = haskey(original_molecules, name) ? casefold(name) : name
        e = RegistryEntry("formula" => name, "nasa" => name, "elements" => elements)
        if !(name in covered)
            get!(reg.molecular, key, e)
            push!(catalogue, CatalogueRow(key, name, "molecular", false, NO_SAFT_REASON, name))
        end
        if data.phase == "S"
            get!(reg.solid, key, e)
            push!(catalogue, CatalogueRow(key, name, "solid", true, "NASA solid", name))
        end
    end
    for (name, e) in reg.solid
        push!(catalogue, CatalogueRow(name, e["formula"], "solid", _has(nasa_db, e["nasa"]), "NASA solid", e["nasa"]))
    end

    # Aqueous ions: charged ePC-SAFT rows joined to SUPCRT records by charge and identity.
    xml = iterxml(databases.supcrt, "Species")
    ready_epc, ready_hkf = ParameterSet[], ParameterSet[]
    charged = [i for i in rows_where(epcsaft, "record_type", "species_parameter")
               if _charge_value(cell(epcsaft, "Charge", i)) != 0]
    groups = OrderedDict{Any,Vector{Int}}()   # groupby('Species', sort=False): NaN names dropped
    for i in charged
        s = cell(epcsaft, "Species", i)
        isna(s) || push!(get!(groups, s, Int[]), i)
    end
    for (raw_name, rows) in groups
        name = String(strip(py_str(raw_name)))
        charge = py_int(first(cell(epcsaft, "Charge", i) for i in rows if !isna(cell(epcsaft, "Charge", i))))
        formula = py_str(first(cell(epcsaft, "Formula", i) for i in rows if !isna(cell(epcsaft, "Formula", i))))
        reg.ion_charge[name] = charge
        reg.ion_saft_formula[name] = formula
        identity = _ion_identity(name, charge)
        matches = [node for node in xml if findtext(node, "Type") == "Aqueous" &&
                   py_float(_truthy(findtext(node, "Charge")) ? findtext(node, "Charge") : 0) == charge &&
                   any(_ion_identity(c, charge) == identity for c in (findtext(node, "Name"), findtext(node, "Formula")))]
        reg.supcrt_ion_aliases[name] = Union{Nothing,String}[findtext(node, "Name") for node in matches]
        epc_data = first_complete_ion_epcsaft_set(epcsaft, name, formula, charge)
        hkf_data = isempty(matches) ? nothing : first_complete_hkf_set(databases.supcrt, name, reg.supcrt_ion_aliases[name])
        k = findfirst(node -> read_hkf_candidate(node) !== nothing, matches)
        elements = k === nothing ? nothing :
                   Elements(m[1] => py_float(m[2]) for m in eachmatch(ELEMENT_PATTERN, something(findtext(matches[k], "Elements"), "")))
        ready = epc_data !== nothing && hkf_data !== nothing && elements !== nothing && !isempty(elements) &&
                all(e -> haskey(AW, e), keys(elements))
        push!(catalogue, CatalogueRow(name, formula, "aqueous ion", ready,
                                      ready ? "ready" : "Complete matching ionic EOS/HKF/element data missing", nothing))
        if ready
            reg.ion[name] = RegistryEntry("formula" => formula, "charge" => charge,
                                          "supcrt" => hkf_data["SUPCRT name"], "elements" => elements)
            push!(ready_epc, epc_data)
            push!(ready_hkf, hkf_data)
        end
    end
    reg.ion_epcsaft_table = ready_epc
    reg.ion_hkf_table = ready_hkf
    # SUPCRT-only records are listed so users can see what exists.
    ion_names = Set(r.Species for r in catalogue if r.Phase == "aqueous ion")
    for node in xml
        name, kind = findtext(node, "Name"), findtext(node, "Type")
        name in ion_names && continue
        push!(catalogue, CatalogueRow(name, findtext(node, "Formula"), "SUPCRT $(py_str(kind))", false,
                                      "SUPCRT reference present; matching EOS/reference adapter required", nothing))
    end
    # drop_duplicates(['Species', 'Phase']), keeping the first
    seen_rows = Set{Tuple{Any,String}}()
    reg.catalogue = [r for r in catalogue if !((r.Species, r.Phase) in seen_rows) && (push!(seen_rows, (r.Species, r.Phase)); true)]
    _register_supcrt_minerals!(reg, databases)
    return reg
end

"""
Single-range Maier-Kelley minerals from supcrt92.xml become pure solids.

Polymorphic records (phase transitions) are listed as unavailable rather
than extrapolated. Mineral <Thermo/HKF> records are not aqueous HKF records.
"""
function _register_supcrt_minerals!(reg::Registry, databases::Databases)
    rows = CatalogueRow[]
    for node in iterxml(databases.supcrt, "Species")
        findtext(node, "Type") == "Mineral" || continue
        name = findtext(node, "Name")
        hkf = findpath(node, "Thermo/HKF")
        reason = ""
        params = OrderedDict{String,Float64}()
        try
            hkf === nothing && throw(ArgumentError("Mineral HKF data missing"))
            ntrans = findtext(hkf, "NumPhaseTrans")
            py_int(_truthy(ntrans) ? ntrans : 0) != 0 &&
                throw(ArgumentError("Polymorphic mineral: transition adapter not implemented"))
            for tag in ("Gf", "Hf", "Sr", "Vr", "a", "b", "c", "Tmax")
                params[tag] = py_float(findtext(hkf, tag))
            end
            if !all(isfinite, values(params)) || params["Vr"] <= 0 || params["Tmax"] < 298.15
                throw(ArgumentError("Invalid/missing mineral parameters"))
            end
            for tag in ("Gf", "Hf", "Sr", "a", "b", "c")
                occursin("cal", getattribute(findchild(hkf, tag), "units", "")) ||
                    throw(ArgumentError("Unexpected units for mineral $tag"))
            end
            getattribute(findchild(hkf, "Vr"), "units") == "cm3/mol" ||
                throw(ArgumentError("Expected Vr in cm3/mol (do not use mislabeled MolarVolume)"))
            elements = Elements(m[1] => py_float(m[2])
                                for m in eachmatch(ELEMENT_PATTERN, something(findtext(node, "Elements"), "")))
            (isempty(elements) || !all(e -> haskey(ATOMIC_WEIGHTS, e), keys(elements))) &&
                throw(ArgumentError("Element composition or atomic weight missing"))
            formula = something(findtext(node, "Formula"), "")
            reg.solid[py_str(name)] = RegistryEntry(
                "formula" => findtext(node, "Formula"), "elements" => elements, "adapter" => "supcrt_mineral",
                "nasa" => nothing, "supcrt_mineral" => params, "source" => "supcrt92.xml",
                "references" => findtext(hkf, "References"), "date" => findtext(hkf, "Date"),
                "hydrate" => occursin('*', formula) && occursin("H2O", formula))
        catch err
            err isa ArgumentError || rethrow()
            reason = err.msg
        end
        push!(rows, CatalogueRow(name, findtext(node, "Formula"), "solid", isempty(reason),
                                 isempty(reason) ? "SUPCRT mineral: single-range Maier-Kelley" : reason, nothing))
    end
    names = Set(r.Species for r in rows)
    kept = [r for r in reg.catalogue if !(r.Species in names && r.Phase in ("SUPCRT Mineral", "solid"))]
    reg.catalogue = vcat(kept, rows)
    return reg
end
