# The model data: databases + species registry (Python: htg/model.py,
# `Model.load` and the species information methods). The calculations
# (`solve`, `sweep`, `check`) come with the equilibrium port.
#
# A Model is built once and then only read, so one Model can be shared by
# all threads.

const REPO_DATA = normpath(joinpath(@__DIR__, "..", "..", "..", "data"))

"""\$HTG_DATA_DIR if set, else the repository's data/ folder."""
default_data_dir() = (d = get(ENV, "HTG_DATA_DIR", ""); isempty(d) ? REPO_DATA : d)

"""
    Model

Thermodynamic model data built from the databases in `data_dir`:
`databases`, `registry`, the `water_model` used for every aqueous property
(SUPCRT92Water() unless chosen otherwise), the user species folder in use
(`user_data`, `nothing`: shipped data only), what that folder contributes
(`user_info`) and how to load the same model again (`load_args`).
"""
struct Model
    databases::Databases
    registry::Registry
    data_dir::Union{Nothing,String}
    water_model::WaterModel
    user_data::Union{Nothing,String}
    user_info::OrderedDict{String,Any}
    load_args::OrderedDict{String,Any}
end

_water_model_name(::SUPCRT92Water) = "SUPCRT92"
_water_model_name(::IAPWS95Water) = "IAPWS95"
_water_model_name(m::WaterModel) = string(typeof(m))

"""
    load_model(data_dir=default_data_dir(); species_file=nothing, nasa_names=nothing,
               reference_temperature_C=350.0, user_data=nothing, water_model=SUPCRT92Water()) -> Model

Read the databases and discover species (Python: `Model.load`).

- `species_file`: registry file (default: data_dir/species.toml)
- `nasa_names`: replaces the [nasa_names] table of the species file
- `reference_temperature_C`: used only to evaluate temperature-dependent SAFT
  parameters while checking which species have complete data
- `user_data`: user species folder (see user_species.jl), read after the
  shipped files; `nothing` (default): shipped data only
- `water_model`: the water model of every aqueous property (one model throughout)
"""
function load_model(data_dir=nothing; species_file=nothing, nasa_names=nothing, reference_temperature_C=350.0,
                    user_data=nothing, water_model::WaterModel=SUPCRT92Water())
    data_dir = data_dir === nothing ? default_data_dir() : String(data_dir)
    databases = load_databases(data_dir)
    config = parsefile_ordered(species_file === nothing ? joinpath(data_dir, "species.toml") : String(species_file))
    nasa_names === nothing || (config["nasa_names"] = OrderedDict{String,Any}(nasa_names))
    user_db, user_info = nothing, OrderedDict{String,Any}()
    if user_data !== nothing
        user_db = load_user_database(user_data)
        user_info = apply_user_layer!(databases, config, user_db)
    end
    registry = Registry(config)
    discover!(registry, databases, Float64(reference_temperature_C) + 273.15)
    user_db === nothing || hide_species!(registry, hidden(user_db))
    user_dir = user_db === nothing ? nothing : user_db.directory
    load_args = OrderedDict{String,Any}(
        "data_dir" => data_dir, "species_file" => species_file === nothing ? nothing : String(species_file),
        "nasa_names" => nasa_names, "reference_temperature_C" => reference_temperature_C,
        "user_data" => user_dir, "water_model" => _water_model_name(water_model))
    return Model(databases, registry, data_dir, water_model, user_dir, user_info, load_args)
end

# ---------------------------------------------------------------- standard states with the model's data
"""Water properties at T [K], P [Pa] from the model's water model."""
water_properties(model::Model, T::Real, P::Real) = water_properties(model.water_model, T, P)

"""g° [J/mol] of aqueous ions at T [K], P [Pa] with the model's HKF records and water model."""
hkf_ion_g0(model::Model, ions, T::Real, P::Real) =
    hkf_ion_g0(model.water_model, model.registry.ion_hkf_table, ions, T, P)

nasa_g0(model::Model, species, T) = nasa_g0(model.registry, model.databases.nasa, species, T)
solid_g0(model::Model, species, T, P) = solid_g0(model.registry, model.databases.nasa, species, T, P)

# ---------------------------------------------------------------- species information
"""
    species(model; include_unavailable=false, phase=nothing, search=nothing) -> Vector{CatalogueRow}

Species catalogue (ready species only, unless `include_unavailable`). `phase`
is a case-insensitive regex on the Phase column; `search` a case-insensitive
substring of Species or Formula.
"""
function species(model::Model; include_unavailable=false, phase=nothing, search=nothing)
    table = available(model.registry; include_unavailable)
    if _truthy(phase)
        pattern = Regex(String(phase), "i")
        table = filter(r -> occursin(pattern, r.Phase), table)
    end
    if _truthy(search)
        text = casefold(String(search))
        hit(x) = x !== nothing && occursin(text, casefold(x))
        table = filter(r -> hit(r.Species) || hit(r.Formula), table)
    end
    return table
end

"""
    mineral_inventory(model; temperature_C=nothing, elements=nothing) -> Vector of rows

SUPCRT minerals available as pure solids, with hydrate flag and T limit
(rows with Species, Formula, Elements, Hydrate, Tmax_K, Within_T_range,
References_in_XML, Date_in_XML).
"""
function mineral_inventory(model::Model; temperature_C=nothing, elements=nothing)
    rows = OrderedDict{String,Any}[]
    for (name, d) in model.registry.solid
        get(d, "adapter", nothing) == "supcrt_mineral" || continue
        if _truthy(elements) && !isempty(elements) && isempty(intersect(Set(elements), keys(d["elements"])))
            continue
        end
        high = d["supcrt_mineral"]["Tmax"]
        push!(rows, OrderedDict{String,Any}(
            "Species" => name, "Formula" => d["formula"], "Elements" => copy(d["elements"]),
            "Hydrate" => d["hydrate"], "Tmax_K" => high,
            "Within_T_range" => temperature_C === nothing || 298.15 <= temperature_C + 273.15 <= high,
            "References_in_XML" => d["references"], "Date_in_XML" => d["date"]))
    end
    return rows
end

"""Python's str() of an exception caught as ValueError/KeyError."""
_exception_text(e::ArgumentError) = e.msg
_exception_text(e::KeyError) = py_repr(e.key isa AbstractString ? e.key : string(e.key))
_exception_text(e) = sprint(showerror, e)

"""
    product_candidates(model, elements; temperature_C=nothing) -> Vector of rows

Neutral gas molecules and pure solids that can carry `elements` (e.g. salt
elements): rows with Species, Phase, Elements, Source, Ready, Reason.
"""
function product_candidates(model::Model, elements; temperature_C=nothing)
    registry, nasa_db = model.registry, model.databases.nasa
    target = Set(elements)
    rows = OrderedDict{String,Any}[]
    for (phase, entries) in (("gas", all_molecular(registry)), ("solid", registry.solid))
        for (name, d) in entries
            (haskey(registry.ion, name) || isempty(intersect(target, keys(d["elements"])))) && continue
            reason = ""
            if get(d, "adapter", nothing) == "supcrt_mineral"
                high = d["supcrt_mineral"]["Tmax"]
                if temperature_C !== nothing && !(298.15 <= temperature_C + 273.15 <= high)
                    reason = "Outside SUPCRT mineral temperature range"
                end
                source = "supcrt92.xml"
            else
                key = get(d, "nasa", nothing)
                record = _has(nasa_db, key) ? nasa_db[key] : nothing
                (record === nothing || record.phase != (phase == "solid" ? "S" : "G")) && continue
                source = "NASA-7 database"
                if temperature_C !== nothing && !(record.T_low <= temperature_C + 273.15 <= record.T_high)
                    reason = "Outside NASA temperature range"
                end
            end
            if phase == "gas"
                T = (temperature_C === nothing ? 350.0 : temperature_C) + 273.15
                try
                    resolve_parameters(registry, model.databases, T, [name])
                catch err
                    (err isa ArgumentError || err isa KeyError) || rethrow()
                    reason = _exception_text(err)
                end
            end
            push!(rows, OrderedDict{String,Any}("Species" => name, "Phase" => phase,
                                                "Elements" => copy(d["elements"]), "Source" => source,
                                                "Ready" => isempty(reason), "Reason" => reason))
        end
    end
    return rows
end
