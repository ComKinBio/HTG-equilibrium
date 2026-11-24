# Case files (Python: tests/test_case_files.py and the case part of
# tests/test_case_format.py; solving is not ported), writing, the ordered
# TOML reader, and - with test/data/cases_python.json - every case file and
# variant against Python field by field.
isdefined(@__MODULE__, :PY_CASES) || include("python_cases.jl")
using HTGEquilibrium: LEGACY_CARBON_KEYS, starts_from_multistart, fluid_phase_count

const CASES_DIR = joinpath(REPO_ROOT, "cases")
const LEGACY_DIR = joinpath(REPO_ROOT, "tests", "data", "legacy_cases")
const TEMPLATE = joinpath(REPO_ROOT, "src", "htg", "templates", "case_template.toml")
const LEGACY_FILES = sort([joinpath(LEGACY_DIR, f) for f in readdir(LEGACY_DIR) if endswith(f, ".toml")])
const CASE_FILES = sort([joinpath(root, f) for (root, _, files) in walkdir(CASES_DIR) for f in files if endswith(f, ".toml")])

new_form(legacy_path) = basename(legacy_path) == "case_template.toml" ? TEMPLATE : joinpath(CASES_DIR, basename(legacy_path))
stem(path) = splitext(basename(path))[1]

"""Everything the calculation takes from the case, as plain values."""
function model_numbers(case)
    feed = case.feed
    dry = dry_element_wt_percent(feed, nothing)
    return Dict("dry" => dry, "water" => total_water_g_per_g_dry(feed), "ash_rest" => 100 - sum(values(dry); init=0.0),
                "salts" => feed.salts, "temperatures" => temperatures(case), "grid" => grid(case),
                "chains" => chains(case), "carbon_bounds" => case.carbon_bounds,
                "solid_carbon_target" => case.solid_carbon_target, "starts" => case.starts)
end

function same_numbers(a, b)
    @test collect(keys(a["dry"])) == collect(keys(b["dry"]))
    for element in keys(a["dry"])
        @test isapprox(a["dry"][element], b["dry"][element]; rtol=1e-12, atol=1e-12)
    end
    @test a["water"] ≈ b["water"] rtol = 1e-12
    @test isapprox(a["ash_rest"], b["ash_rest"]; rtol=1e-12, atol=1e-12)
    for key in ("salts", "temperatures", "grid", "chains", "carbon_bounds", "solid_carbon_target", "starts")
        @test a[key] == b[key]
    end
end

"""The case without the feed and temperature list (compared by number elsewhere)."""
rest_of_case(case) = replace(case; feed=Feed(), feed_v2=nothing, temperature_C=nothing, refine=nothing,
                             comparison=nothing, description="")

legacy(; changes...) = begin
    data = read_toml(joinpath(LEGACY_DIR, "case_template.toml"))
    for (section, values) in changes
        values === nothing ? delete!(data, string(section)) : (data[string(section)] = values)
    end
    data
end
D(pairs...) = OrderedDict{String,Any}(pairs...)
message(f) = try f(); "" catch e; e isa HTGEquilibrium.ValueError ? e.msg : rethrow() end

@testset "case files (htg/case.py)" begin
    @testset "TOML key order" begin
        text = """
        b = 1
        a = { z = 1, y = { x = 2, w = 3 } }
        "quoted key" = 'not = a key'
        dotted.second = 1
        dotted.first = 2
        list = [ { q = 1, p = 2 },
                 # a comment with key = value
                 { o = 3 } ]
        text = \"\"\"
        [not.a.table]
        k = v\"\"\"

        [t.inner]   # comment [x]
        c = 1
        b = 2

        [t]
        z = 0
        y = [1, 2.5, "s", [3]]

        [[arr]]
        n = 1
        m = 2
        [[arr]]
        l = 3
        [arr.sub]
        k = 4
        """
        d = parse_toml(text)
        @test collect(keys(d)) == ["b", "a", "quoted key", "dotted", "list", "text", "t", "arr"]
        @test collect(keys(d["a"])) == ["z", "y"] && collect(keys(d["a"]["y"])) == ["x", "w"]
        @test collect(keys(d["dotted"])) == ["second", "first"]
        @test collect(keys(d["list"][1])) == ["q", "p"]
        @test d["text"] == "[not.a.table]\nk = v" && d["quoted key"] == "not = a key"
        @test collect(keys(d["t"])) == ["inner", "z", "y"] && collect(keys(d["t"]["inner"])) == ["c", "b"]
        @test d["t"]["y"] == Any[1, 2.5, "s", Any[3]] && d["t"]["y"][1] isa Int64
        @test collect(keys(d["arr"][1])) == ["n", "m"] && collect(keys(d["arr"][2])) == ["l", "sub"]
        @test d isa OrderedDict{String,Any}
        @test_throws Exception parse_toml("a = ")
    end

    @testset "all case files load" begin
        @test !isempty(CASE_FILES)
        for path in CASE_FILES
            case = load_case(path)
            @test case.name == stem(path)
            if py_truth(case.comparison)
                @test isfile(case.comparison["file"])
            end
            @test is_v2(read_toml(path))
        end
        @test load_case(TEMPLATE).name == "my_case"
    end

    @testset "old and converted files give the same case: $(stem(path))" for path in LEGACY_FILES
        old = load_case(path)
        shipped = load_case(new_form(path))
        converted = case_from_dict(legacy_to_v2(read_toml(path)); name=stem(path))
        for new in (shipped, converted)
            @test new.feed_v2 !== nothing
            same_numbers(model_numbers(old), model_numbers(new))
            @test rest_of_case(new) == rest_of_case(old)
        end
        if old.feed.ash_wt_percent !== nothing  # ash given in the old file is kept
            @test converted.feed.ash_wt_percent == old.feed.ash_wt_percent
        end
    end

    @testset "legacy carbon and starts" begin
        data = legacy(constraints=D("gas_carbon_min_fraction" => 0.9, "solid_carbon_fraction" => 0.05,
                                    "phase_carbon_fractions" => D("organic_liquid" => 0.01)),
                      solver=D("method" => "ipopt", "multistart" => 4))
        v2 = legacy_to_v2(data)
        @test v2["format"] == 2 && first(keys(v2)) == "format"
        @test v2["constraints"] == D("carbon" => D("gas" => D("min" => 0.9), "solid" => D("exact" => 0.05),
                                                   "organic_liquid" => D("exact" => 0.01)))
        @test collect(keys(v2["constraints"]["carbon"])) == ["gas", "solid", "organic_liquid"]
        # gas, aqueous and organic liquid have species (the cut-off at 390 °C is above 300 °C)
        @test v2["solver"] == D("method" => "ipopt", "starts" => 1 + 2 * 3 + 4)
        @test legacy_to_v2(legacy(solver=D("multistart" => 0)))["solver"] == D("starts" => 1)
        hot = legacy(solver=D("multistart" => 2))
        hot["conditions"]["temperature_C"] = [400, 450]  # aqueous removed at every temperature
        @test legacy_to_v2(hot)["solver"] == D("starts" => 1 + 2 * 2 + 2)
        @test case_from_dict(data).starts == 11 && case_from_dict(hot).starts == 7
        @test legacy_to_v2(v2) == v2  # already format 2
    end

    @testset "refine is merged into the temperatures" begin
        data = legacy(sweep=D("continuation" => true, "refine" => D("start_C" => 340, "stop_C" => 370, "step_C" => 10)))
        v2 = legacy_to_v2(data)
        @test v2["conditions"]["temperature_C"] == [300, 340, 350, 360, 370, 400, 450]
        @test all(t -> t isa Int, v2["conditions"]["temperature_C"])
        @test v2["sweep"] == D("continuation" => true)
        @test temperatures(case_from_dict(v2)) == temperatures(case_from_dict(data))
        only = legacy_to_v2(legacy(sweep=D("refine" => D("start_C" => 340, "stop_C" => 370, "step_C" => 10))))
        @test !haskey(only, "sweep")
    end

    @testset "legacy_carbon_bounds mapping" begin
        @test legacy_carbon_bounds(gas_carbon_fraction=0.8) == (D("gas" => D("exact" => 0.8)), "char")
        @test legacy_carbon_bounds(gas_carbon_min_fraction=0.9, phase_carbon_fractions=D("aqueous" => 0.0)) ==
              (D("gas" => D("min" => 0.9), "aqueous" => D("exact" => 0.0)), "char")
        @test legacy_carbon_bounds(phase_carbon_fractions=D("solid" => 0.1)) == (D("solid" => D("exact" => 0.1)), "all_solids")
        @test occursin("not both", message(() -> legacy_carbon_bounds(gas_carbon_fraction=0.8, gas_carbon_min_fraction=0.9)))
        @test occursin("Conflicting", message(() -> legacy_carbon_bounds(gas_carbon_fraction=0.8,
                                                                       phase_carbon_fractions=D("gas" => 0.7))))
    end

    @testset "the case moves old settings into the new ones" begin
        case = replace(load_case(TEMPLATE); gas_carbon_fraction=0.8, phase_carbon_fractions=D("aqueous" => 0.0),
                       multistart=2)
        @test case.carbon_bounds == D("gas" => D("exact" => 0.8), "aqueous" => D("exact" => 0.0))
        @test case.gas_carbon_fraction === nothing && case.phase_carbon_fractions === nothing
        @test case.starts == 1 + 2 * 3 + 2 && case.multistart === nothing
        @test replace(case) == case
        @test occursin("given twice", message(() -> replace(case; gas_carbon_min_fraction=0.5)))
        @test_throws ArgumentError replace(case; no_such_field=1)
    end

    @testset "carbon bound validation: $msg" for (bounds, msg) in [
        (D("gas" => D("exact" => 0.5, "min" => 0.2)), "Gas carbon: exact cannot be combined with min or max."),
        (D("solid" => D("min" => 0.4, "max" => 0.2)), "Solid: the minimum is larger than the maximum."),
        (D("aqueous" => D("max" => 1.5)), "Aqueous carbon: max must be a fraction between 0 and 1"),
        (D("organic_liquid" => D("at_least" => 0.1)), "Organic liquid carbon: use exact, or min and/or max"),
        (D("liquid" => D("exact" => 0.1)), "Unknown phase 'liquid'")]
        @test occursin(msg, message(() -> check_carbon_bounds(bounds)))
        @test occursin(msg, message(() -> validate(replace(load_case(TEMPLATE); carbon_bounds=bounds))))
    end

    @testset "carbon bounds in a file" begin
        dir = mktempdir()
        text = replace(read(TEMPLATE, String), "# [constraints.carbon]\n# gas = { min = 0.9 }\n# solid = { max = 0.05 }" =>
                       "[constraints.carbon]\ngas = { min = 0.9 }\nsolid = { max = 0.05 }")
        write(joinpath(dir, "c.toml"), text)
        case = load_case(joinpath(dir, "c.toml"))
        @test case.carbon_bounds == D("gas" => D("min" => 0.9), "solid" => D("max" => 0.05))
        write(joinpath(dir, "bad.toml"), replace(text, "solid = { max = 0.05 }" => "solid = { min = 0.2, max = 0.05 }"))
        @test message(() -> load_case(joinpath(dir, "bad.toml"))) == "Solid: the minimum is larger than the maximum."
    end

    @testset "starts validation" begin
        case = load_case(TEMPLATE)
        @test case.starts == 1 && case.distinct_composition_tol == 1e-4 && case.distinct_gibbs_tol == 1.0
        for bad in (0, 2.5, "x")
            @test occursin("starts", message(() -> validate(replace(case; starts=bad))))
        end
        @test validate(replace(case; starts=3.0)).starts === 3
        @test occursin("either starts or the old multistart", message(() -> replace(case; starts=3, multistart=1)))
    end

    @testset "starts_from_multistart" begin
        @test starts_from_multistart(0, 3) == 1
        @test starts_from_multistart(4, 2) == 9
        @test_throws CaseError("cannot convert float NaN to integer") starts_from_multistart(NaN, 2)
        @test_throws CaseError("cannot convert float infinity to integer") starts_from_multistart(Inf, 2)
        phases = D("gas" => ["water"], "aqueous" => ["water"], "organic_liquid" => String[], "solid" => ["carbon"])
        @test fluid_phase_count(phases) == 2
        @test fluid_phase_count(phases, false) == 1
        @test fluid_phase_count(phases, true, 390, [300, 400]) == 2
        @test fluid_phase_count(phases, true, 390, [390, 400]) == 1
    end

    @testset "water model ([aqueous] water_model, Julia only)" begin
        data = read_toml(TEMPLATE)
        @test !haskey(data["aqueous"], "water_model")
        default = case_from_dict(data)
        @test default.water_model == "SUPCRT92" && water_model_of(default) isa SUPCRT92Water
        for (name, type) in ("SUPCRT92" => SUPCRT92Water, "IAPWS-95" => IAPWS95Water)
            data["aqueous"]["water_model"] = name
            case = case_from_dict(data)
            @test case.water_model == name && water_model_of(case) isa type
        end
        data["aqueous"]["water_model"] = "IAPWS95"
        want = "water_model must be one of ['SUPCRT92', 'IAPWS-95'], not 'IAPWS95'."
        @test message(() -> case_from_dict(data)) == want
        @test message(() -> water_model_of("IAPWS95")) == want
        @test message(() -> validate(replace(default; water_model=nothing))) ==
              "water_model must be one of ['SUPCRT92', 'IAPWS-95'], not None."
        # written only when not the default
        @test !occursin("water_model", to_toml(default))
        iapws = replace(default; water_model="IAPWS-95")
        @test occursin("\n[aqueous]\n", to_toml(iapws)) && occursin("\nwater_model = \"IAPWS-95\"\n", to_toml(iapws))
        @test case_from_dict(parse_toml(to_toml(iapws))) == iapws
    end

    @testset "writing: load, save, load gives the same case: $(stem(path))" for path in [CASE_FILES; TEMPLATE]
        case = load_case(path)
        written = joinpath(mktempdir(), basename(path))
        HTGEquilibrium.save(case, written)
        @test read(written, String) == to_toml(case) == dump_case(case)
        @test load_case(written) == case
    end

    @testset "writing old files: format 2, same model numbers: $(stem(path))" for path in LEGACY_FILES
        case = load_case(path)
        again = case_from_dict(parse_toml(to_toml(case)); name=case.name)
        @test parse_toml(to_toml(case))["format"] == 2
        same_numbers(model_numbers(again), model_numbers(case))
    end

    @testset "writing: a fresh file, comments are not kept" begin
        commented(text) = any(l -> startswith(lstrip(l), "#") || occursin(r"\s#", l), split(text, '\n'))
        original = read(TEMPLATE, String)
        @test commented(original)                         # the template is full of comments ...
        written = to_toml(load_case(TEMPLATE))
        @test !commented(written)                         # ... the written file has none
        @test !any(commented(to_toml(load_case(p))) for p in CASE_FILES)
        @test startswith(written, "format = 2\nname = \"my_case\"\n")
    end

    @testset "mixed feed keys are refused" begin
        dir = mktempdir()
        write(joinpath(dir, "c.toml"), replace(read(TEMPLATE, String), "total_kg = 5.0" => "total_kg = 5.0\nwater_g_per_g_dry = 4.0"))
        @test occursin("mixes the new keys", message(() -> load_case(joinpath(dir, "c.toml"))))
        data = read_toml(TEMPLATE)
        data["feed"] = D("elements_wt_percent" => data["feed"]["elements_wt_percent"])
        @test occursin("total_kg and solids", message(() -> case_from_dict(data)))
        delete!(data, "format")  # an old file: elements only, no water
        @test case_from_dict(data).feed.water_g_per_g_dry == 0
    end

    @testset "species feed in a file" begin
        text = read(TEMPLATE, String)
        start, stop = findfirst("total_kg = 5.0", text)[1], findfirst("# Or, instead", text)[1]
        text = text[1:prevind(text, start)] * text[stop:end]
        text = replace(text, "# [feed.species]\n# water = { value = 4.5, unit = \"kg\" }      # always present; \"kg\" or \"mol\"\n" *
                             "# cellulose = { value = 0.5, unit = \"kg\" }" =>
                             "[feed.species]\nwater = { value = 4.5, unit = \"kg\" }\ncellulose = { value = 0.5, unit = \"kg\" }")
        text = replace(text, "[feed]\n# The wet feed going into the reactor: solids plus water.\n" => "")
        dir = mktempdir()
        write(joinpath(dir, "c.toml"), text)
        case = load_case(joinpath(dir, "c.toml"))
        @test case.feed_v2 == D("species" => D("water" => D("value" => 4.5, "unit" => "kg"),
                                               "cellulose" => D("value" => 0.5, "unit" => "kg")))
        elements = species_elements_function(D("water" => D("H" => 2, "O" => 1), "cellulose" => D("C" => 6, "H" => 10, "O" => 5)))
        mm = molar_mass_function(elements)
        r = resolved(case.feed, mm)
        @test r.water_g_per_g_dry ≈ 9.0
        @test r.species_moles["cellulose"] ≈ 500 / mm("cellulose")
    end

    @testset "grid, chains, at" begin
        case = load_case(TEMPLATE)
        @test grid(case) == [(T, 240.0) for T in [300.0, 350.0, 390.0, 400.0, 450.0]]
        @test temperatures(case) == [300.0, 350.0, 390.0, 400.0, 450.0]   # the cut-off 390 °C lies inside
        @test chains(case) == [grid(case)]
        multi = replace(case; pressure_bar=[250, 240.0])
        @test pressures(multi) == [240.0, 250.0]
        @test chains(multi) == [[(T, P) for T in temperatures(case)] for P in (240.0, 250.0)]
        @test length(chains(replace(multi; continuation=false))) == 10
        @test occursin("several pressures", message(() -> at(multi, 300)))
        point = at(case, 375)
        @test point.temperature_C === 375.0 && point.pressure_bar === 240.0 && is_single_point(point)
        pts = replace(case; temperature_C=nothing, pressure_bar=nothing, points=[[300, 240], [250, 250], [300, 240]])
        @test grid(pts) == [(300.0, 240.0), (250.0, 250.0)] && chains(pts) == [grid(pts)]
        @test temperatures(pts) == [250.0, 300.0]
        @test effective_carbon_bounds(D("aqueous" => D("min" => 0.1), "gas" => D("max" => 0.5)), true) ==
              (D("gas" => D("max" => 0.5)), ["aqueous carbon"])
        @test BOUND_KINDS == ("exact", "min", "max") && PHASE_LABELS["organic_liquid"] == "Organic liquid"
    end

    if PY_CASES !== nothing
        @testset "against Python: $(length(PY_CASES["files"])) case files" begin
            for want in PY_CASES["files"]
                @test_same file_dump(want["path"]) want want["path"]
            end
        end
        @testset "against Python: $(length(PY_CASES["variants"])) variants" begin
            for want in PY_CASES["variants"]
                text = variant_text(want["base"], [(old, new) for (old, new) in want["replace"]])
                got = text_dump(text, "variant")
                @test_same got strip_keys(want, "label", "base", "replace") want["label"]
            end
        end
        @info "Case/feed floats compared with Python (rtol $FLOAT_RTOL, atol $FLOAT_ATOL): $(FLOATS.compared); " *
              "$(FLOATS.inexact) not identical; largest relative deviation $(FLOATS.worst), " *
              "largest absolute deviation of the rest $(FLOATS.worst_abs)"
    end
end
