# The user species database, as far as loading it into the model goes
# (Python: htg/user_species.py: resolve_user_dir, UserDatabase.load, owners,
# hidden, apply_user_layer, hide_species). Editing, importing and the GUI
# helpers (save_species, SpeciesSpec, check_species, ...) come later (M4).
#
# Folder (default platformdirs.user_data_dir("htg-equilibrium")/species;
# override with $HTG_USER_DATA, where the folder holds the files directly):
#
#     user_pcsaft.csv     PC-SAFT rows, same columns as data/PROJECT_PCSAFT.csv
#     user_epcsaft.csv    ePC-SAFT rows, same columns as data/PROJECT_EPCSAFT.csv
#     user_thermo.txt     NASA-7 records (CHEMKIN 4-line format)
#     user_supcrt.json    aqueous HKF and mineral Maier-Kelley records (fields of supcrt92.xml)
#     user_species.toml   [species.<name>] index, [nasa_names], [solid.<name>],
#                         [[binary_interaction]], [removed]
#
# `load_model(...; user_data=DIR)` reads data/ first, then these files
# (`apply_user_layer!`). A user record with the same species name replaces the
# shipped one; the shipped files are never written.

const APP_NAME = "htg-equilibrium"
const USER_FILES = OrderedDict("pcsaft" => "user_pcsaft.csv", "epcsaft" => "user_epcsaft.csv",
                               "nasa" => "user_thermo.txt", "supcrt" => "user_supcrt.json",
                               "species" => "user_species.toml")

const PCSAFT_COLUMNS = [
    "record_type", "Species", "Formula", "Type", "Charge", "MW_g_mol", "Model", "m", "sigma_A",
    "u_over_kB_K", "Association_scheme", "epsilonAB_over_kB_K", "kappaAB", "Dipole_D",
    "Dipole_fraction", "Quadrupole_DA", "Quadrupole_fraction", "Ionization_energy_eV",
    "Source_extra_value", "Species2", "kij", "Temperature dependence", "T_unit", "K1", "K2",
    "RefID", "Authors", "Title", "Journal", "Year", "DOI", "Notes", "Source_table"]
const EPCSAFT_COLUMNS = [
    "record_type", "Species", "Formula", "Type", "Charge", "Species2", "Model", "Parameter",
    "Value", "Units", "kij", "Temperature dependence", "lij", "kij_hb",
    "Dielectric correlation/value", "Variable/type", "RefID", "Authors", "Title", "Journal",
    "Year", "DOI", "Notes"]

# SUPCRT92 fields and units (as in supcrt92.xml; calories).
const HKF_FIELDS = ["Gf", "Hf", "Sr", "a1", "a2", "a3", "a4", "c1", "c2", "wref"]
const MAIER_KELLEY_FIELDS = ["Gf", "Hf", "Sr", "Vr", "a", "b", "c", "Tmax"]
const SUPCRT_UNITS = Dict("Gf" => "cal/mol", "Hf" => "cal/mol", "Sr" => "cal/(mol*K)", "a1" => "cal/(mol*bar)",
                          "a2" => "cal/mol", "a3" => "(cal*K)/(mol*bar)", "a4" => "(cal*K)/mol",
                          "c1" => "cal/(mol*K)", "c2" => "(cal*K)/mol", "wref" => "cal/mol", "Vr" => "cm3/mol",
                          "Tmax" => "K", "a" => "cal/(mol*K)", "b" => "cal/(mol*K^2)", "c" => "(cal*K)/mol")

# ---------------------------------------------------------------- folder
"""
    default_user_dir()

\$HTG_USER_DATA if set, else the platform's user data folder for
'htg-equilibrium' (as Python's platformdirs) plus /species.
"""
function default_user_dir()
    env = get(ENV, "HTG_USER_DATA", "")
    isempty(env) || return expanduser(env)
    base = if Sys.iswindows()
        joinpath(get(ENV, "LOCALAPPDATA", joinpath(homedir(), "AppData", "Local")), APP_NAME, APP_NAME)
    elseif Sys.isapple()
        joinpath(homedir(), "Library", "Application Support", APP_NAME)
    else
        xdg = get(ENV, "XDG_DATA_HOME", "")
        joinpath(isempty(strip(xdg)) ? joinpath(homedir(), ".local", "share") : xdg, APP_NAME)
    end
    return joinpath(base, "species")
end

"""The user database folder for a command: `nothing` when switched off (--no-user-species)."""
function resolve_user_dir(path=nothing; enabled=true)
    enabled || return nothing
    return (path === nothing || path == "") ? default_user_dir() : expanduser(String(path))
end

# ---------------------------------------------------------------- the database
"""
    UserDatabase

The user's species files in one folder: `pcsaft`, `epcsaft` (DataTables of
text), `nasa` (record name => 4-line record text), `supcrt` (user_supcrt.json
records), `config` (user_species.toml).
"""
mutable struct UserDatabase
    directory::String
    pcsaft::DataTable
    epcsaft::DataTable
    nasa::OrderedDict{String,String}
    supcrt::Vector{Any}
    config::OrderedDict{String,Any}
end

function _normalize_config!(c)
    get!(c, "species", OrderedDict{String,Any}())
    get!(c, "nasa_names", OrderedDict{String,Any}())
    get!(c, "solid", OrderedDict{String,Any}())
    get!(c, "binary_interaction", Any[])
    get!(get!(c, "removed", OrderedDict{String,Any}()), "species", Any[])
    return c
end

UserDatabase(directory) = UserDatabase(String(directory), empty_table(PCSAFT_COLUMNS), empty_table(EPCSAFT_COLUMNS),
                                       OrderedDict{String,String}(), Any[],
                                       _normalize_config!(OrderedDict{String,Any}()))

user_path(db::UserDatabase, key) = joinpath(db.directory, USER_FILES[key])

"""
    load_user_database(directory) -> UserDatabase

Read the folder; missing files count as empty (the folder need not exist).
Python: `UserDatabase.load`.
"""
load_user_database(db::UserDatabase) = db
function load_user_database(directory)
    db = UserDatabase(directory)
    isfile(user_path(db, "pcsaft")) && (db.pcsaft = read_text_csv(user_path(db, "pcsaft"), PCSAFT_COLUMNS))
    isfile(user_path(db, "epcsaft")) && (db.epcsaft = read_text_csv(user_path(db, "epcsaft"), EPCSAFT_COLUMNS))
    isfile(user_path(db, "nasa")) && (db.nasa = _split_nasa_text(read_text(user_path(db, "nasa"))))
    if isfile(user_path(db, "supcrt"))
        text = read(user_path(db, "supcrt"), String)
        data = isempty(strip(text)) ? Dict() : JSON.parse(text)
        db.supcrt = collect(Any, get(data, "species", Any[]))
    end
    if isfile(user_path(db, "species"))
        db.config = _normalize_config!(read_toml(user_path(db, "species")))
    end
    return db
end

"""{record name => record text} of a CHEMKIN thermo file (strict columns)."""
function _split_nasa_text(text)
    lines = py_splitlines(text; keepends=true)
    out = OrderedDict{String,String}()
    i = 1
    while i + 3 <= length(lines)
        group = [rstrip(line, ('\r', '\n')) for line in lines[i:i+3]]
        if all(length(line) >= 80 && py_slice(line, 79, 80) == string(k) for (k, line) in enumerate(group))
            out[String(strip(py_slice(group[1], 0, 18)))] = join(group, "\n") * "\n"
            i += 4
        else
            i += 1
        end
    end
    return out
end

_molecule_name(value) = casefold(strip(py_str(value)))
_user_key(name) = normalize_ion_name(name)

function _number(value)
    v = py_tryfloat(value)
    return (v === nothing || !isfinite(v)) ? nothing : v
end
_finite(value) = (v = py_tryfloat(value); v !== nothing && isfinite(v))

"""{species name => set of user files with records for it}."""
function owners(db::UserDatabase)
    out = OrderedDict{String,Set{String}}()
    function add(name, key)
        (name === nothing || isna(name) || py_str(name) == "") && return
        push!(get!(out, py_str(name), Set{String}()), USER_FILES[key])
    end
    record_owner = Dict(get(entry, "nasa_record", nothing) => name for (name, entry) in db.config["species"]
                        if _truthy(get(entry, "nasa_record", nothing)))
    for i in 1:nrows(db.pcsaft)
        py_str(cell(db.pcsaft, "record_type", i)) == "pure_component" &&
            add(_molecule_name(cell(db.pcsaft, "Species", i)), "pcsaft")
    end
    for i in 1:nrows(db.epcsaft)
        if py_str(cell(db.epcsaft, "record_type", i)) == "species_parameter"
            charged = !(_number(cell(db.epcsaft, "Charge", i)) in (nothing, 0.0))
            species = cell(db.epcsaft, "Species", i)
            add(charged ? String(strip(py_str(species))) : _molecule_name(species), "epcsaft")
        end
    end
    for record in keys(db.nasa)
        add(get(record_owner, record, record), "nasa")
    end
    for record in db.supcrt
        add(get(record, "Name", nothing), "supcrt")
    end
    for name in vcat(collect(keys(db.config["species"])), collect(keys(db.config["nasa_names"])),
                     collect(keys(db.config["solid"])))
        add(name, "species")
    end
    return out
end

species_names(db::UserDatabase) = sort(collect(keys(owners(db))); by=casefold)
hidden(db::UserDatabase) = collect(db.config["removed"]["species"])

# ---------------------------------------------------------------- SUPCRT records
_elements_text(elements) = join("$e($(_g(py_float(n))))" for (e, n) in elements)

"""XML <Species> element (as in supcrt92.xml) from a user_supcrt.json record."""
function supcrt_node(record)
    node = XMLNode("Species")
    function sub(parent, tag, text, units=nothing)
        child = XMLNode(tag; text=py_str(text))
        units === nothing || (child.attributes["units"] = units)
        push!(parent.children, child)
        return child
    end
    sub(node, "Name", record["Name"])
    formula = get(record, "Formula", nothing)
    sub(node, "Formula", _truthy(formula) ? formula : record["Name"])
    elements = get(record, "Elements", nothing)
    elements = (elements === nothing || isempty(elements)) ? Dict() : elements
    sub(node, "Elements", elements isa AbstractString ? elements : _elements_text(elements))
    kind = get(record, "Type", "Aqueous")
    sub(node, "Type", kind)
    charge = get(record, "Charge", nothing)
    kind == "Aqueous" && sub(node, "Charge", py_int(_truthy(charge) && charge != 0 ? charge : 0))
    thermo = XMLNode("Thermo")
    push!(node.children, thermo)
    hkf = XMLNode("HKF")
    push!(thermo.children, hkf)
    sub(hkf, "References", get(record, "References", "user"))
    sub(hkf, "Date", get(record, "Date", ""))
    values = get(record, "HKF", nothing)
    values = values === nothing ? Dict() : values
    for name in (kind == "Mineral" ? MAIER_KELLEY_FIELDS : HKF_FIELDS)
        value = get(values, name, nothing)
        (value === nothing || !_finite(value)) && continue
        sub(hkf, name, py_float(value), SUPCRT_UNITS[name])
    end
    if kind == "Mineral"
        n = get(values, "NumPhaseTrans", nothing)
        sub(hkf, "NumPhaseTrans", py_int(_truthy(n) && n != 0 ? n : 0))
    end
    return node
end

# ---------------------------------------------------------------- merging into the model
"""
    apply_user_layer!(databases, config, user_db) -> info

Merge the user records into `databases` and the species.toml `config` (in
place). Returns {species name => {"files" => [...], "overrides" => Bool}}: what
the user database contributes, and whether it replaces shipped records.
"""
function apply_user_layer!(databases::Databases, config::AbstractDict, user_db::UserDatabase)
    owner_files = owners(user_db)
    shipped_pc, shipped_epc = databases.pcsaft, databases.epcsaft
    shipped_nasa = Set(keys(databases.nasa))
    shipped_supcrt = Set(findtext(node, "Name") for node in iterxml(databases.supcrt, "Species"))
    shipped_pc_names = Set(_user_key(cell(shipped_pc, "Species", i)) for i in rows_where(shipped_pc, "record_type", "pure_component"))
    shipped_epc_names = Set(_user_key(cell(shipped_epc, "Species", i)) for i in rows_where(shipped_epc, "record_type", "species_parameter"))
    shipped_solids = Set(keys(get(config, "solid", Dict())))

    # SAFT tables: user rows replace shipped rows of the same species and come first,
    # so they win where the lookup goes by formula.
    pc_user = rows_where(user_db.pcsaft, "record_type", "pure_component")
    epc_user = rows_where(user_db.epcsaft, "record_type", "species_parameter")
    names = union(Set(_user_key(cell(user_db.pcsaft, "Species", i)) for i in pc_user),
                  Set(_user_key(cell(user_db.epcsaft, "Species", i)) for i in epc_user))
    if !isempty(pc_user)
        databases.pcsaft = _merge_frame(shipped_pc, user_db.pcsaft, pc_user, "pure_component", names,
                                        USER_FILES["pcsaft"])
    end
    if !isempty(epc_user) || !isempty(pc_user)
        databases.epcsaft = _merge_frame(shipped_epc, user_db.epcsaft, epc_user, "species_parameter", names,
                                         USER_FILES["epcsaft"])
    end

    # NASA-7: same record name replaces.
    for text in values(user_db.nasa)
        merge!(databases.nasa, read_nasa_lines(py_splitlines(text; keepends=true)))
    end

    # SUPCRT92: user records replace same names and go first.
    if !isempty(user_db.supcrt)
        root = databases.supcrt
        user_names = Set(get(r, "Name", nothing) for r in user_db.supcrt)
        filter!(n -> !(n.tag == "Species" && findtext(n, "Name") in user_names), root.children)
        position = something(findfirst(n -> n.tag == "Species", root.children), length(root.children) + 1)
        for (offset, record) in enumerate(user_db.supcrt)
            insert!(root.children, position + offset - 1, supcrt_node(record))
        end
    end

    # species.toml tables.
    nasa_names = get!(config, "nasa_names", OrderedDict{String,Any}())
    for (name, record) in user_db.config["nasa_names"]
        _has(databases.nasa, record) && (nasa_names[name] = record)
    end
    solids = get!(config, "solid", OrderedDict{String,Any}())
    for (name, entry) in user_db.config["solid"]
        all(k -> haskey(entry, k), ("formula", "nasa", "elements")) && (solids[name] = copy(entry))
    end
    pairs = get!(config, "binary_interaction", Any[])
    for item in user_db.config["binary_interaction"]
        components = something(get(item, "components", nothing), [])
        (length(components) != 2 || !_finite(get(item, "kij", nothing))) && continue
        pair = Set(components)
        filter!(p -> Set(p["components"]) != pair, pairs)
        push!(pairs, OrderedDict{String,Any}("components" => collect(components), "kij" => py_float(item["kij"])))
    end

    info = OrderedDict{String,Any}()
    for (name, files) in owner_files
        key = _user_key(name)
        overrides = key in shipped_pc_names || key in shipped_epc_names || name in shipped_nasa ||
                    name in shipped_supcrt || name in shipped_solids
        record = get(get(user_db.config["species"], name, Dict()), "nasa_record", nothing)
        overrides = overrides || record in shipped_nasa
        info[name] = OrderedDict{String,Any}("files" => sort(collect(files)), "overrides" => overrides)
    end
    return info
end

"""
User rows first (index -n..-1, with a 'source_file' column), then the shipped
rows not replaced by them; the shipped numeric columns stay numeric
(`pd.to_numeric(errors='coerce')`). User columns not in the shipped table are dropped.
"""
function _merge_frame(shipped::DataTable, user::DataTable, user_rows, record_type, names, source_file)
    kept = [i for i in 1:nrows(shipped)
            if !(py_str(cell(shipped, "record_type", i)) == record_type &&
                 _user_key(cell(shipped, "Species", i)) in names)]
    columns = unique(vcat(shipped.columns, "source_file"))
    data = Dict{String,Vector{Any}}()
    for c in columns
        user_values = c == "source_file" ? Any[source_file for _ in user_rows] : Any[cell(user, c, i) for i in user_rows]
        shipped_values = Any[cell(shipped, c, i) for i in kept]
        if c in shipped.numeric
            user_values = Any[(v = isna(v) ? nothing : py_tryfloat(v); v === nothing ? NaN : v) for v in user_values]
        end
        data[c] = vcat(user_values, shipped_values)
    end
    index = vcat(collect(-length(user_rows):-1), shipped.index[kept])
    return DataTable(columns, data, copy(shipped.numeric), index)
end

"""Drop hidden shipped species from a discovered registry and its catalogue."""
function hide_species!(reg::Registry, names)
    names = Set(names)
    isempty(names) && return reg
    for table in (reg.molecular, reg.organic, reg.solid, reg.ion)
        for name in intersect(names, keys(table))
            delete!(table, name)
        end
    end
    filter!(r -> !(r.Species in names), reg.catalogue)
    return reg
end
