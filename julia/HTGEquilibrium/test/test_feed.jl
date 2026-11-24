# Feed arithmetic of case-file format 2 and the conversion of old feeds
# (Python: the feed part of tests/test_case_format.py); with
# test/data/cases_python.json also feed_v2 functions on many [feed] tables
# against Python.
isdefined(@__MODULE__, :PY_CASES) || include("python_cases.jl")
using HTGEquilibrium: _plus_one, refined_temperatures, UNMODELLED_MESSAGE

const D_ = OrderedDict{String,Any}
# Elements of a few species, for species-amount feeds without the databases.
const FEW_ELEMENTS = D_("water" => D_("H" => 2, "O" => 1), "cellulose" => D_("C" => 6, "H" => 10, "O" => 5),
                        "methanol" => D_("C" => 1, "H" => 4, "O" => 1))
const few_species = species_elements_function(FEW_ELEMENTS)
const FEED = D_("total_kg" => 10.0, "solids" => D_("value" => 15.0, "unit" => "%"),
                "ash" => D_("value" => 8.45, "unit" => "%"),
                "elements_wt_percent" => D_("C" => 30.4, "H" => 4.41, "O" => 53.58, "N" => 2.66, "S" => 0.51))
with(d, pairs...) = (out = deepcopy(d); for (k, v) in pairs; out[k] = v; end; out)
feed_error(f) = try f(); nothing catch e; e isa FeedError ? e.msg : rethrow() end
r2(x, n) = round(x; digits=n)

@testset "feed (htg/feed_v2.py)" begin
    @testset "legacy feed rules" begin
        dry = legacy_feed_to_v2(D_("elements_wt_percent" => D_("C" => 45.0, "H" => 5.0, "O" => 40.0),
                                   "ash_wt_percent" => 10.0, "water_g_per_g_dry" => 3.0, "moisture_wt_percent" => 20.0))
        @test dry["total_kg"] ≈ 1 + 3.0 + 20 / 80
        @test dry["solids"] == D_("value" => 1.0, "unit" => "kg")
        @test dry["ash"] == D_("value" => 10.0, "unit" => "%")
        @test all(isapprox.(collect(values(dry["elements_wt_percent"])), [50.0, 50 / 9, 400 / 9]))
        daf = legacy_feed_to_v2(D_("elements_wt_percent" => D_("C" => 50.0, "H" => 6.0, "O" => 44.0),
                                   "basis" => "dry_ash_free", "ash_wt_percent" => 5.0))
        @test daf["elements_wt_percent"] == D_("C" => 50.0, "H" => 6.0, "O" => 44.0)
        @test daf["ash"] == D_("value" => 5.0, "unit" => "%") && daf["total_kg"] == 1.0
        @test !haskey(legacy_feed_to_v2(D_("elements_wt_percent" => D_("C" => 50.0, "H" => 6.0, "O" => 44.0))), "ash")
        err = try
            legacy_feed_to_v2(D_("elements_wt_percent" => D_("C" => 50.0, "H" => 6.0, "O" => 34.0),
                                 "ash_wt_percent" => 5.0, "unmodelled_wt_percent" => 5.0))
        catch e
            e
        end
        @test err isa LegacyConversionError && err isa HTGEquilibrium.ValueError && occursin("unmodelled_wt_percent", err.msg)
        err = try legacy_feed_to_v2(D_("species_moles" => D_("cellulose" => 1.0), "water_g_per_g_dry" => 4.0)) catch e; e end
        @test err isa LegacyConversionError && occursin("species data", err.msg)
    end

    @testset "legacy species feed" begin
        mm = molar_mass_function(few_species)
        v2 = legacy_feed_to_v2(D_("species_moles" => D_("cellulose" => 2.0, "methanol" => 1.0), "water_g_per_g_dry" => 4.0), mm)
        dry_g = 2 * mm("cellulose") + mm("methanol")
        @test collect(keys(v2["species"])) == ["water", "cellulose", "methanol"]
        @test v2["species"]["water"]["value"] ≈ 4.0 * dry_g / 1000 && v2["species"]["water"]["unit"] == "kg"
        @test v2["species"]["cellulose"] == D_("value" => 2.0, "unit" => "mol")
        @test v2["species"]["methanol"] == D_("value" => 1.0, "unit" => "mol")
        model = feed_v2_to_model(v2; molar_mass=mm)
        @test model.species_moles == D_("cellulose" => 2.0, "methanol" => 1.0)
        @test model.water_g_per_g_dry ≈ 4.0 rtol = 1e-14
        deferred = feed_v2_to_model(v2)
        @test deferred.species_amounts == v2["species"]
        err = try total_water_g_per_g_dry(deferred) catch e; e end
        @test err isa CaseError && occursin("resolved", err.msg)
        @test resolved(deferred, mm) == model
    end

    @testset "feed summary matches the Case tab example" begin
        summary = feed_summary(FEED)
        @test summary["mode"] == "elements"
        @test D_(e => r2(v, 2) for (e, v) in summary["used"]) == D_("C" => 33.20, "H" => 4.82, "O" => 58.52, "N" => 2.91, "S" => 0.56)
        @test r2(summary["entered_sum"], 2) == 91.56
        @test r2(summary["water_kg"], 3) == 8.500 && r2(summary["water_mol"], 1) == 471.8
        @test r2(summary["ash_kg"], 3) == 0.127 && r2(summary["ash_free_kg"], 3) == 1.373
        @test r2(summary["elements"]["C"]["kg"], 3) == 0.456
        @test r2(summary["elements"]["C"]["mol"], 2) == 37.96
        @test [row["name"] for row in summary["rows"]] == ["Total", "Water", "Solids", "Ash", "Ash-free solids",
                                                         "C", "H", "O", "N", "S"]
        @test summary["rows"][1] == D_("name" => "Total", "kg" => 10.0, "mol" => nothing)
        kg = feed_summary(with(FEED, "solids" => D_("value" => 1.5, "unit" => "kg"), "ash" => D_("value" => 0.12675, "unit" => "kg")))
        @test kg["ash_free_kg"] ≈ summary["ash_free_kg"]
    end

    @testset "feed summary, species mode" begin
        feed = D_("species" => D_("water" => D_("value" => 4.5, "unit" => "kg"), "cellulose" => D_("value" => 0.5, "unit" => "kg"),
                                  "methanol" => D_("value" => 2.0, "unit" => "mol")))
        summary = feed_summary(feed; species_elements=few_species)
        M = molar_mass_function(few_species)
        @test summary["mode"] == "species"
        names = [row["name"] for row in summary["rows"]]
        @test names[1:5] == ["Total", "Water", "Solids (not water)", "cellulose", "methanol"]
        @test Set(names[6:end]) == Set(["C", "H", "O"])
        dry = 0.5 + 2.0 * M("methanol") / 1000
        @test summary["solids_kg"] ≈ dry
        @test summary["total_kg"] ≈ 4.5 + dry
        @test summary["elements"]["C"]["mol"] ≈ 6 * 500 / M("cellulose") + 2.0
        @test occursin("species data", feed_error(() -> feed_summary(feed)))
    end

    @testset "feed errors: $msg" for (change, msg) in [
        ("ash" => D_("value" => 2.0, "unit" => "kg"), "Ash can't be more than the solids."),
        ("elements_wt_percent" => D_(), "Elements: give at least one element."),
        ("solids" => D_("value" => 15.0, "unit" => "g"), "Solids: the unit must be"),
        ("total_kg" => 0, "Total mass must be positive."),
        ("solids" => D_("value" => 12.0, "unit" => "kg"), "Solids can't be more than the total mass.")]
        @test startswith(feed_error(() -> feed_summary(with(FEED, change))), msg)
        @test startswith(feed_error(() -> feed_v2_to_model(with(FEED, change))), msg)
    end

    @testset "feed_v2_to_model" begin
        model = feed_v2_to_model(FEED)
        @test model.basis == "dry_ash_free" && model.ash_wt_percent == 8.45 && model.moisture_wt_percent == 0
        @test model.water_g_per_g_dry ≈ 8.5 / 1.5
        @test sum(values(model.elements_wt_percent)) ≈ 100
        @test model.elements_wt_percent["C"] ≈ 30.4 / 91.56 * 100
        @test sum(values(dry_element_wt_percent(model, nothing))) ≈ 100 - 8.45
        @test feed_v2_to_model(with(FEED, "ash" => D_("value" => 0.12675, "unit" => "kg"))).ash_wt_percent ≈ 8.45
        @test feed_v2_to_model(D_(k => v for (k, v) in FEED if k != "ash")).ash_wt_percent === nothing
        s, a = feed_fractions(model)
        @test s ≈ 0.15 && a ≈ 0.0845
        salts = feed_v2_to_model(with(FEED, "salts" => D_("KCl" => D_("mol_per_g_dry_feed" => 1e-5, "elements" => D_("K" => 1, "Cl" => 1)))))
        @test salts.salts["KCl"] isa Salt && salts.salts["KCl"].mol_per_g_dry_feed == 1e-5
        @test elemental_masses(FEED)["ash_free"] ≈ 1.5 * (1 - 0.0845)
        @test ash_percent(FEED) == 8.45 && ash_percent(with(FEED, "ash" => nothing)) === nothing
    end

    @testset "feed fractions of old feeds" begin
        old = load_case(joinpath(REPO_ROOT, "tests", "data", "legacy_cases", "yakaboylu2013_fig7_8.toml")).feed
        s, a = feed_fractions(old)
        @test s ≈ 1 / 4.999935304 && a ≈ 0.0844688764
    end

    @testset "_plus_one" begin
        @test _plus_one(3.999935304) == 4.999935304
        @test _plus_one(0.1) == 1.1 && _plus_one(1e-30) == 1.0 && _plus_one(-0.0) == 1.0
    end

    if PY_CASES !== nothing
        @testset "against Python: feed_v2 on $(length(PY_CASES["feeds"])) [feed] tables" begin
            mm = case_molar_mass()
            for want in PY_CASES["feeds"]
                feed = decode(want["feed"])
                got = OrderedDict{String,Any}(
                    "feed" => enc(feed),
                    "feed_summary" => enc(attempt(feed_summary, deepcopy(feed), nothing, SPECIES_ELEMENTS)),
                    "feed_summary_no_species" => enc(attempt(feed_summary, deepcopy(feed))),
                    "feed_v2_to_model" => enc(attempt(feed_v2_to_model, deepcopy(feed))),
                    "feed_v2_to_model_mm" => enc(attempt(feed_v2_to_model, deepcopy(feed); molar_mass=mm)),
                    "elemental_masses" => enc(attempt(elemental_masses, deepcopy(feed))),
                    "ash_percent" => enc(attempt(ash_percent, deepcopy(feed))))
                if haskey(got["feed_v2_to_model"], "ok")
                    got["feed_fractions"] = enc(attempt(feed_fractions, feed_v2_to_model(deepcopy(feed)), mm))
                end
                @test_same got want repr(feed)
            end
        end
        @testset "against Python: legacy_feed_to_v2 on $(length(PY_CASES["legacy_feeds"])) old feeds" begin
            mm = case_molar_mass()
            for want in PY_CASES["legacy_feeds"]
                feed = decode(want["feed"])
                got = OrderedDict{String,Any}("feed" => enc(feed),
                                              "legacy_feed_to_v2" => enc(attempt(legacy_feed_to_v2, deepcopy(feed))),
                                              "legacy_feed_to_v2_mm" => enc(attempt(legacy_feed_to_v2, deepcopy(feed), mm)))
                @test_same got want repr(feed)
            end
        end
        @testset "against Python: _plus_one, refined_temperatures" begin
            @test all(_plus_one(x) === Float64(y) for (x, y) in PY_CASES["plus_one"])
            for (g, refine, want) in PY_CASES["refined_temperatures"]
                @test_same enc(refined_temperatures(decode(g), decode(refine))) want "refine $refine"
            end
        end
    end
end
