# App preferences (Python: htg/settings.py): display choices only, never
# anything that changes a calculation.
#
# Stored as TOML in `settings.toml` in the config folder, which is, in order:
# the `config_dir` argument of `load_settings`, the `HTG_CONFIG_DIR`
# environment variable, or the folder `platformdirs.user_config_dir(
# "htg-equilibrium", appauthor=False)` gives in Python (see `user_config_dir`).
# Example file:
#
#     temperature_unit = "C"        # Units.TEMPERATURE_UNITS
#     pressure_unit = "bar"         # Units.PRESSURE_UNITS
#     color_theme = "Default"
#     line_width = 2
#     input_formats_open = false
#
#     [export_defaults.one]         # Save plot dialog, one entry per layout
#     width_mm = 85
#     height_mm = 65
#     dpi = 300
#     format = "png"                # png, svg or pdf
#     shared_legend = true
#
# Layouts: `one` (one panel), `grid` (all four, 2 x 2), `row` (all four in a
# row), `present` (phases present, in a row).
#
#     settings = load_settings()                     # once, at start-up (Python: Settings.load())
#     settings["pressure_unit"]                      # "bar"         (settings.get('pressure_unit'))
#     settings["pressure_unit"] = "MPa"              # validates and writes the file (settings.set(...))
#     export_defaults(settings, "row")               # OrderedDict("width_mm" => 250, ...) (a copy)
#     set_export_defaults(settings, "row"; width_mm=240)
#
# A missing, unreadable or invalid file (or a single invalid value) falls back
# to the defaults with a logged warning. Failing to write is also only a
# warning. Writing keeps the key order of the defaults and drops comments (as
# Python's tomlkit.dumps of the values); the file is the same text Python
# writes for the same values.

const SETTINGS_APP_NAME = "htg-equilibrium"
const SETTINGS_FILE_NAME = "settings.toml"
const SETTINGS_ENV_VAR = "HTG_CONFIG_DIR"

const EXPORT_LAYOUTS = ("one", "grid", "row", "present")
const EXPORT_FORMATS = ("png", "svg", "pdf")

_export_default(width_mm, height_mm) = OrderedDict{String,Any}(
    "width_mm" => width_mm, "height_mm" => height_mm, "dpi" => 300, "format" => "png", "shared_legend" => true)

"""The default preferences (Python: `settings.DEFAULTS`); copy before changing."""
const SETTINGS_DEFAULTS = OrderedDict{String,Any}(
    "temperature_unit" => "C",
    "pressure_unit" => "bar",
    "color_theme" => "Default",
    "line_width" => 2,
    "input_formats_open" => false,
    "export_defaults" => OrderedDict{String,Any}(
        "one" => _export_default(85, 65), "grid" => _export_default(170, 130),
        "row" => _export_default(250, 65), "present" => _export_default(190, 65)),
)

# Python's isinstance(v, (int, float)) and not isinstance(v, bool): TOML gives Int64 or Float64.
_positive_number(v) = (v isa Integer || v isa AbstractFloat) && !(v isa Bool) && v > 0
_positive_int(v) = v isa Integer && !(v isa Bool) && v > 0

const _SETTING_CHECKS = OrderedDict{String,Function}(
    "temperature_unit" => v -> v isa AbstractString && v in Units.TEMPERATURE_UNITS,
    "pressure_unit" => v -> v isa AbstractString && v in Units.PRESSURE_UNITS,
    "color_theme" => v -> v isa AbstractString && v != "",
    "line_width" => _positive_number,
    "input_formats_open" => v -> v isa Bool,
)
const _EXPORT_CHECKS = OrderedDict{String,Function}(
    "width_mm" => _positive_number,
    "height_mm" => _positive_number,
    "dpi" => _positive_int,
    "format" => v -> v isa AbstractString && v in EXPORT_FORMATS,
    "shared_legend" => v -> v isa Bool,
)

"""An unknown setting or export layout (Python: KeyError)."""
struct SettingsKeyError <: Exception
    msg::String
end
Base.showerror(io::IO, e::SettingsKeyError) = print(io, e.msg)

"""An invalid value for a setting (Python: ValueError)."""
struct SettingsValueError <: ValueError
    msg::String
end
Base.showerror(io::IO, e::SettingsValueError) = print(io, e.msg)

"""
    user_config_dir() -> String

`platformdirs.user_config_dir("htg-equilibrium", appauthor=False)` (platformdirs
4.14): Linux `\$XDG_CONFIG_HOME/htg-equilibrium` when that is an absolute path,
else `~/.config/htg-equilibrium`; macOS the same `\$XDG_CONFIG_HOME` rule, else
`~/Library/Application Support/htg-equilibrium`; Windows
`%LOCALAPPDATA%\\htg-equilibrium` (no author folder, because appauthor=False;
`WIN_PD_OVERRIDE_LOCAL_APPDATA` first, as platformdirs; Julia reads the
LOCALAPPDATA variable where platformdirs asks Windows for the same folder).
"""
function user_config_dir()
    if Sys.iswindows()
        override = strip(get(ENV, "WIN_PD_OVERRIDE_LOCAL_APPDATA", ""))
        base = if _fully_qualified_windows(override)
            override
        else
            local_appdata = get(ENV, "LOCALAPPDATA", "")
            isempty(local_appdata) ? joinpath(homedir(), "AppData", "Local") : local_appdata
        end
        return joinpath(normpath(base), SETTINGS_APP_NAME)
    end
    xdg = get(ENV, "XDG_CONFIG_HOME", "")
    startswith(xdg, "/") && return joinpath(xdg, SETTINGS_APP_NAME)   # posixpath.isabs
    base = Sys.isapple() ? joinpath(homedir(), "Library", "Application Support") : joinpath(homedir(), ".config")
    return joinpath(base, SETTINGS_APP_NAME)
end

# platformdirs' _is_fully_qualified: a drive letter with a separator, or a UNC path.
_fully_qualified_windows(path) =
    occursin(r"^[A-Za-z]:[\\/]", path) || occursin(r"^[\\/]{2}[^\\/]+[\\/][^\\/]+", path)

"""Path of settings.toml: in `config_dir`, `\$HTG_CONFIG_DIR` or the user config folder."""
function settings_path(config_dir=nothing)
    config_dir === nothing || return joinpath(String(config_dir), SETTINGS_FILE_NAME)
    env = get(ENV, SETTINGS_ENV_VAR, "")
    isempty(env) || return joinpath(env, SETTINGS_FILE_NAME)
    return joinpath(user_config_dir(), SETTINGS_FILE_NAME)
end

"""
    Settings(path, values=copy of SETTINGS_DEFAULTS)

Preferences read from and written to one settings file (Python: `Settings`).
Read with `load_settings`; `settings[key]`, `settings[key] = value`,
`export_defaults`, `set_export_defaults`, `as_dict`, `save`.
"""
mutable struct Settings
    path::String
    values::OrderedDict{String,Any}
end
Settings(path) = Settings(String(path), deepcopy(SETTINGS_DEFAULTS))

Base.show(io::IO, s::Settings) = print(io, "Settings(", repr(s.path), ")")

"""
    load_settings(config_dir=nothing) -> Settings

Read `settings.toml` from the config folder; defaults for anything missing or
invalid (Python: `Settings.load`). Each ignored value is a logged warning.
"""
function load_settings(config_dir=nothing)
    path = settings_path(config_dir)
    values = deepcopy(SETTINGS_DEFAULTS)
    ispath(path) || return Settings(path, values)
    stored = try
        text = read(path, String)
        isvalid(text) || error("the file is not UTF-8 text")
        parse_toml(text)
    catch exc
        exc isa InterruptException && rethrow()
        reason = exc isa TOML.ParserError ? sprint(showerror, exc) : (exc isa ErrorException ? exc.msg : sprint(showerror, exc))
        @warn "Could not read $path ($(_one_line(reason))); using the default settings."
        return Settings(path, values)
    end
    for (key, value) in stored
        if key == "export_defaults"
            _merge_export!(values["export_defaults"], value, path)
        elseif !haskey(_SETTING_CHECKS, key)
            @warn "Ignoring unknown setting \"$key\" in $path."
        elseif _SETTING_CHECKS[key](value)
            values[key] = value
        else
            @warn "Ignoring invalid value $(py_repr(value)) for \"$key\" in $path."
        end
    end
    return Settings(path, values)
end

_one_line(text) = join(split(strip(text), '\n'; keepempty=false), " ")

function _merge_export!(target, stored, path)
    if !(stored isa AbstractDict)
        @warn "Ignoring invalid \"export_defaults\" in $path."
        return nothing
    end
    for (layout, values) in stored
        if !haskey(target, layout) || !(values isa AbstractDict)
            @warn "Ignoring export_defaults.$layout in $path."
            continue
        end
        for (key, value) in values
            if haskey(_EXPORT_CHECKS, key) && _EXPORT_CHECKS[key](value)
                target[layout][key] = value
            else
                @warn "Ignoring export_defaults.$layout.$key = $(py_repr(value)) in $path."
            end
        end
    end
    return nothing
end

"""A preference value (a copy; export_defaults: use `export_defaults(settings, layout)`)."""
function Base.getindex(s::Settings, key::AbstractString)
    haskey(SETTINGS_DEFAULTS, key) || throw(SettingsKeyError("Unknown setting \"$key\""))
    return deepcopy(s.values[key])
end
Base.get(s::Settings, key::AbstractString) = s[key]

"""Change a preference and write the file if the value changed (Python: `settings.set(key, value)`)."""
function Base.setindex!(s::Settings, value, key::AbstractString)
    haskey(_SETTING_CHECKS, key) ||
        throw(SettingsKeyError("Unknown setting \"$key\"" * (key == "export_defaults" ? " (use set_export_defaults)" : "")))
    _SETTING_CHECKS[key](value) || throw(SettingsValueError("Invalid value $(py_repr(value)) for setting \"$key\""))
    if !_same_value(s.values[key], value)
        s.values[key] = value isa AbstractString ? String(value) : value
        save(s)
    end
    return s
end

# Python's != : 2 == 2.0, True is never stored where a number is (the checks refuse Bool).
_same_value(a, b) = a == b && (a isa Bool) == (b isa Bool)

"""`(width_mm, height_mm, dpi, format, shared_legend)` for a layout, as an OrderedDict (a copy)."""
function export_defaults(s::Settings, layout::AbstractString)
    layout in EXPORT_LAYOUTS ||
        throw(SettingsKeyError("Unknown export layout \"$layout\"; use one of $(join(EXPORT_LAYOUTS, ", "))"))
    return copy(s.values["export_defaults"][layout])
end

"""Remember Save plot values for a layout, e.g. `set_export_defaults(settings, "one"; width_mm=90)`."""
function set_export_defaults(s::Settings, layout::AbstractString; values...)
    current = export_defaults(s, layout)
    for (key, value) in values
        k = String(key)
        haskey(_EXPORT_CHECKS, k) || throw(SettingsKeyError("Unknown export setting \"$k\""))
        _EXPORT_CHECKS[k](value) || throw(SettingsValueError("Invalid value $(py_repr(value)) for export setting \"$k\""))
    end
    if any(!_same_value(current[String(k)], v) for (k, v) in values)
        for (k, v) in values
            s.values["export_defaults"][layout][String(k)] = v isa AbstractString ? String(v) : v
        end
        save(s)
    end
    return s
end

"""All values (a copy)."""
as_dict(s::Settings) = deepcopy(s.values)

"""Write the file (all values). Errors are logged, not raised."""
function save(s::Settings)
    try
        mkpath(dirname(s.path))
        temporary = splitext(s.path)[1] * ".toml.tmp"
        write(temporary, settings_toml(s.values))
        mv(temporary, s.path; force=true)
    catch exc
        exc isa InterruptException && rethrow()
        exc isa Base.IOError || exc isa SystemError || exc isa ArgumentError || rethrow()
        @warn "Could not save the settings to $(s.path) ($(_one_line(sprint(showerror, exc))))."
    end
    return nothing
end

# ---------------------------------------------------------------- writing (tomlkit.dumps)
"""The settings file text for `values`, as Python's `tomlkit.dumps(values)` writes it."""
function settings_toml(values::AbstractDict)
    io = IOBuffer()
    tables = Pair{String,Any}[]
    for (key, value) in values
        if value isa AbstractDict
            push!(tables, key => value)
        else
            println(io, _tomlkit_key(key), " = ", _tomlkit_value(value))
        end
    end
    for (name, table) in tables
        _tomlkit_table(io, name, table)
    end
    return String(take!(io))
end

function _tomlkit_table(io, name, table)
    scalars = [(k, v) for (k, v) in table if !(v isa AbstractDict)]
    if !isempty(scalars)
        print(io, "\n[", name, "]\n")
        for (k, v) in scalars
            println(io, _tomlkit_key(k), " = ", _tomlkit_value(v))
        end
    end
    for (k, v) in table
        v isa AbstractDict && _tomlkit_table(io, name * "." * _tomlkit_key(k), v)
    end
end

_tomlkit_key(key) = occursin(r"^[A-Za-z0-9_-]+$", key) ? key : _tomlkit_string(key)
_tomlkit_value(v::Bool) = v ? "true" : "false"
_tomlkit_value(v::Integer) = string(v)
_tomlkit_value(v::AbstractFloat) = isnan(v) ? "nan" : isinf(v) ? (v > 0 ? "inf" : "-inf") : py_repr(v)
_tomlkit_value(v::AbstractString) = _tomlkit_string(v)

const _TOMLKIT_ESCAPES = Dict('\b' => "\\b", '\t' => "\\t", '\n' => "\\n", '\f' => "\\f", '\r' => "\\r",
                           '"' => "\\\"", '\\' => "\\\\")
function _tomlkit_string(s)
    io = IOBuffer()
    print(io, '"')
    for c in s
        if haskey(_TOMLKIT_ESCAPES, c)
            print(io, _TOMLKIT_ESCAPES[c])
        elseif c < ' ' || c == '\x7f'
            print(io, "\\u", lpad(string(UInt32(c); base=16), 4, '0'))
        else
            print(io, c)
        end
    end
    print(io, '"')
    return String(take!(io))
end
