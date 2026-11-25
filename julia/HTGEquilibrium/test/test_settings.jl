# App preferences (Python: tests/test_settings.py), and files written and read
# by the Python app (test/data/python_settings, written by test/dump_settings.py).
using HTGEquilibrium, Test, JSON
using HTGEquilibrium: OrderedDict
using HTGEquilibrium: as_dict, settings_toml, user_config_dir

const SETTINGS_DATA = joinpath(@__DIR__, "data", "python_settings")

"""Run f() and return (its value, the messages of the warnings it logged)."""
function with_warnings(f)
    logs, value = Test.collect_test_logs(f; min_level=Base.CoreLogging.Warn)
    return value, [string(l.message) for l in logs]
end

"""JSON values as the settings hold them (ordered, Int vs Float kept)."""
plain(x::AbstractDict) = OrderedDict{String,Any}(String(k) => plain(v) for (k, v) in x)
plain(x) = x

@testset "settings (htg/settings.py)" begin
    @testset "defaults when missing" begin
        dir = mktempdir()
        prefs = load_settings(dir)
        @test prefs.path == joinpath(dir, "settings.toml")
        @test as_dict(prefs) == SETTINGS_DEFAULTS
        @test prefs["temperature_unit"] == "C" && prefs["pressure_unit"] == "bar"
        @test prefs["color_theme"] == "Default" && prefs["line_width"] === 2
        @test prefs["input_formats_open"] === false
        @test export_defaults(prefs, "one") == OrderedDict("width_mm" => 85, "height_mm" => 65, "dpi" => 300,
                                                           "format" => "png", "shared_legend" => true)
        sizes = Dict(layout => (export_defaults(prefs, layout)["width_mm"], export_defaults(prefs, layout)["height_mm"])
                     for layout in EXPORT_LAYOUTS)
        @test sizes == Dict("one" => (85, 65), "grid" => (170, 130), "row" => (250, 65), "present" => (190, 65))
        @test !ispath(prefs.path)   # nothing written until a value changes
        @test collect(keys(SETTINGS_DEFAULTS)) ==
              ["temperature_unit", "pressure_unit", "color_theme", "line_width", "input_formats_open", "export_defaults"]
    end

    @testset "values survive a restart" begin
        dir = mktempdir()
        prefs = load_settings(dir)
        prefs["pressure_unit"] = "MPa"
        prefs["temperature_unit"] = "K"
        prefs["line_width"] = 1.5
        prefs["input_formats_open"] = true
        set_export_defaults(prefs, "row"; width_mm=240, format="pdf")
        again = load_settings(dir)
        @test again["pressure_unit"] == "MPa" && again["temperature_unit"] == "K"
        @test again["line_width"] == 1.5 && again["input_formats_open"] === true
        @test export_defaults(again, "row") == OrderedDict("width_mm" => 240, "height_mm" => 65, "dpi" => 300,
                                                           "format" => "pdf", "shared_legend" => true)
        @test export_defaults(again, "one") == SETTINGS_DEFAULTS["export_defaults"]["one"]
        stored = HTGEquilibrium.TOML.parsefile(again.path)
        @test stored["export_defaults"]["row"]["width_mm"] == 240
        @test get(again, "pressure_unit") == "MPa"   # Python's settings.get(key)
    end

    @testset "environment variable" begin
        dir = mktempdir()
        withenv("HTG_CONFIG_DIR" => joinpath(dir, "env")) do
            @test settings_path() == joinpath(dir, "env", "settings.toml")
            @test settings_path(dir) == joinpath(dir, "settings.toml")
            load_settings()["color_theme"] = "Set2"
        end
        @test load_settings(joinpath(dir, "env"))["color_theme"] == "Set2"
    end

    @testset "user config folder (platformdirs)" begin
        withenv("HTG_CONFIG_DIR" => nothing, "XDG_CONFIG_HOME" => nothing) do
            if Sys.iswindows()
                @test settings_path() == joinpath(ENV["LOCALAPPDATA"], "htg-equilibrium", "settings.toml")
            elseif Sys.isapple()
                @test settings_path() == joinpath(homedir(), "Library", "Application Support", "htg-equilibrium",
                                                  "settings.toml")
            else
                @test settings_path() == joinpath(homedir(), ".config", "htg-equilibrium", "settings.toml")
            end
        end
        if !Sys.iswindows()
            withenv("XDG_CONFIG_HOME" => "/somewhere/config") do
                @test user_config_dir() == "/somewhere/config/htg-equilibrium"
            end
            withenv("XDG_CONFIG_HOME" => "relative/config") do   # not absolute: ignored, as platformdirs
                @test !startswith(user_config_dir(), "relative")
            end
        end
    end

    @testset "writes only on change" begin
        prefs = load_settings(mktempdir())
        prefs["pressure_unit"] = "bar"
        set_export_defaults(prefs, "one"; width_mm=85)
        set_export_defaults(prefs, "one"; width_mm=85.0)   # Python: 85 == 85.0
        prefs["line_width"] = 2.0
        @test !ispath(prefs.path)
    end

    @testset "invalid values are refused" begin
        prefs = load_settings(mktempdir())
        @test_throws SettingsValueError prefs["pressure_unit"] = "psi"
        @test_throws SettingsValueError prefs["line_width"] = 0
        @test_throws SettingsValueError prefs["line_width"] = true
        @test_throws SettingsValueError prefs["input_formats_open"] = 1
        @test_throws SettingsKeyError prefs["solver"] = "ipopt"
        @test_throws SettingsKeyError prefs["export_defaults"] = Dict()
        @test_throws SettingsKeyError prefs["solver"]
        @test_throws SettingsKeyError export_defaults(prefs, "poster")
        @test_throws SettingsValueError set_export_defaults(prefs, "one"; format="jpg")
        @test_throws SettingsValueError set_export_defaults(prefs, "one"; dpi=300.0)
        @test_throws SettingsKeyError set_export_defaults(prefs, "one"; colour="red")
        @test !ispath(prefs.path)
        @test sprint(showerror, SettingsValueError("x")) == "x"
        @test SettingsValueError("x") isa ValueError
        err = try prefs["export_defaults"] = 1 catch e e end
        @test sprint(showerror, err) == "Unknown setting \"export_defaults\" (use set_export_defaults)"
        err = try prefs["line_width"] = -1 catch e e end
        @test sprint(showerror, err) == "Invalid value -1 for setting \"line_width\""
    end

    @testset "returned values are copies" begin
        prefs = load_settings(mktempdir())
        export_defaults(prefs, "one")["width_mm"] = 1
        prefs["export_defaults"]["one"]["width_mm"] = 1
        @test export_defaults(prefs, "one")["width_mm"] == 85
    end

    @testset "broken file gives defaults" begin
        dir = mktempdir()
        write(joinpath(dir, "settings.toml"), "pressure_unit = \"MPa\n[[[")
        prefs, warnings = with_warnings(() -> load_settings(dir))
        @test as_dict(prefs) == SETTINGS_DEFAULTS
        @test length(warnings) == 1 && occursin("Could not read", warnings[1])
        write(joinpath(dir, "settings.toml"), [0xff, 0xfe, 0x00])   # not UTF-8 (Python: UnicodeDecodeError)
        prefs, warnings = with_warnings(() -> load_settings(dir))
        @test as_dict(prefs) == SETTINGS_DEFAULTS && occursin("Could not read", only(warnings))
    end

    @testset "bad values are ignored one by one" begin
        dir = mktempdir()
        write(joinpath(dir, "settings.toml"),
              "pressure_unit = \"MPa\"\ntemperature_unit = \"F\"\nline_width = \"thick\"\nunknown = 1\n" *
              "export_defaults = { one = { width_mm = 90, dpi = -1 }, poster = { width_mm = 1 } }\n")
        prefs, warnings = with_warnings(() -> load_settings(dir))
        @test prefs["pressure_unit"] == "MPa"
        @test prefs["temperature_unit"] == "C" && prefs["line_width"] == 2
        @test export_defaults(prefs, "one")["width_mm"] == 90 && export_defaults(prefs, "one")["dpi"] == 300
        text = join(warnings, "\n")
        for word in ["temperature_unit", "line_width", "unknown", "dpi", "poster"]
            @test occursin(word, text)
        end
    end

    @testset "unwritable folder only warns" begin
        dir = mktempdir()
        blocker = joinpath(dir, "file")
        write(blocker, "")
        prefs = load_settings(joinpath(blocker, "config"))
        _, warnings = with_warnings(() -> (prefs["pressure_unit"] = "atm"))
        @test prefs["pressure_unit"] == "atm"
        @test length(warnings) == 1 && occursin("Could not save", warnings[1])
    end

    @testset "file written by Python" begin
        expected = plain(JSON.parsefile(joinpath(SETTINGS_DATA, "written.json"); dicttype=OrderedDict))
        python_text = read(joinpath(SETTINGS_DATA, "written", "settings.toml"), String)
        prefs, warnings = with_warnings(() -> load_settings(joinpath(SETTINGS_DATA, "written")))
        @test isempty(warnings)
        @test as_dict(prefs) == expected
        @test prefs["line_width"] isa Float64 && export_defaults(prefs, "row")["width_mm"] isa Int
        @test prefs["color_theme"] == "Set2 \"bright\" \\ é\ttab"
        # Julia writes the same text for the same values, so Python reads Julia's file as its own.
        @test settings_toml(as_dict(prefs)) == python_text
        copy_dir = mktempdir()
        copied = load_settings(copy_dir)
        for key in ("temperature_unit", "pressure_unit", "color_theme", "line_width", "input_formats_open")
            copied[key] = prefs[key]
        end
        for layout in EXPORT_LAYOUTS
            set_export_defaults(copied, layout; (Symbol(k) => v for (k, v) in export_defaults(prefs, layout))...)
        end
        @test read(copied.path, String) == python_text
    end

    @testset "hand-written file read as Python reads it" begin
        expected = JSON.parsefile(joinpath(SETTINGS_DATA, "mixed.json"); dicttype=OrderedDict)
        path = joinpath(SETTINGS_DATA, "mixed", "settings.toml")
        prefs, warnings = with_warnings(() -> load_settings(dirname(path)))
        @test as_dict(prefs) == plain(expected["values"])
        @test [replace(w, path => "<path>") for w in warnings] == expected["warnings"]
        # writing drops the comments and uses the default key order
        dir = mktempdir()
        cp(path, joinpath(dir, "settings.toml"))
        mixed = with_warnings(() -> load_settings(dir))[1]
        mixed["line_width"] = 3
        text = read(mixed.path, String)
        @test !occursin('#', text) && startswith(text, "temperature_unit = \"C\"\npressure_unit = \"MPa\"\n")
        @test occursin("line_width = 3\n", text)
    end

    @testset "TOML strings and numbers as tomlkit writes them" begin
        S = HTGEquilibrium
        @test S._tomlkit_string("a\tb\nc") == "\"a\\tb\\nc\""
        @test S._tomlkit_string("x\x01y\x7fz") == "\"x\\u0001y\\u007fz\""
        @test S._tomlkit_string("é ✓ it's") == "\"é ✓ it's\""
        @test [S._tomlkit_value(v) for v in (1e16, 1e-5, 0.1, 2.0, Inf, 123456789012345678, 1.5e300)] ==
              ["1e+16", "1e-05", "0.1", "2.0", "inf", "123456789012345678", "1.5e+300"]
        # every written form reads back to the same value
        for v in (1e16, 1e-5, 0.1, 2.0, Inf, 1.5e300, "a\tb\"\\\x01é")
            @test HTGEquilibrium.TOML.parse("k = " * S._tomlkit_value(v))["k"] == v
        end
    end
end
