# The Case tab without the page (Python: tests/test_case_tab.py, test_case_format.py
# for the documents): the wording and numbers of case_form.jl, CaseDocument, and
# the CaseEditor driven as the page drives it (edit!, share_unit!, ...). The page
# itself is checked by the self-test (test/selftest.jl).
using HTGApp: check_feed, conflict_notice, element_error, element_symbol, entered_sum_notice, file_status, format_kg,
              format_mol, initial_guesses_tooltip, parse_number, point_count, point_count_text, problem_text,
              band_points, band_hits, band_notice, FieldError, LEGACY_NOTICE, SALTS_NOTICE, LOADING_MESSAGE,
              case_from_data, document_case, save_document, dirty, has_salts, notices, CaseEditor, form_view, edit!,
              share_unit!, add_element!, remove_element!, carbon_choice!, set_unit!, switch_conditions!,
              set_input_formats!, switch_feed!, water_unit!, add_species!, species_name!, remove_species!, open!,
              save!, write!, new!, queue!, check_case, check_done!, set_choice!, set_water_model!, status_line,
              first_problem, with_extension, case_tab_html, WATER_MODEL_MOLECULAR, WATER_MODEL_NO_AQUEOUS
using HTGEquilibrium: preflight_messages

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const FEED = OrderedDict{String,Any}(
    "total_kg" => 10.0, "solids" => OrderedDict{String,Any}("value" => 15.0, "unit" => "%"),
    "ash" => OrderedDict{String,Any}("value" => 8.45, "unit" => "%"),
    "elements_wt_percent" => OrderedDict{String,Any}("C" => 30.4, "H" => 4.41, "O" => 53.58, "N" => 2.66, "S" => 0.51))

"""A stand-in for supcrt92_band_check: `hits` (T °C) are in a band at every pressure."""
fake_bands(hits; calls=Ref(0)) = function (T, P; margin_C=0.5)
    calls[] += 1
    Ps = P isa AbstractVector ? P : fill(P, length(T))
    return [(T_C=t, P_bar=p, status=t in hits ? :in_band : :ok, band=nothing) for (t, p) in zip(T, Ps)]
end

"""An editor over a new document (or `document`), the band check synchronous."""
function editor(state; document=nothing, band_check=fake_bands(Float64[]))
    state.document = something(document, new_document("form_test"))
    return CaseEditor(state; async_bands=false, band_check)
end

texts(ed) = ed.texts
errors(ed) = form_view(ed)["errors"]

@testset "Case tab" begin
    @testset "feed numbers (WP4 acceptance)" begin
        summary, field, message = check_feed(FEED)
        @test message === nothing
        @test [HTGApp.format_percent(summary["used"][e]) for e in ("C", "H", "O", "N", "S")] == ["33.20", "4.82", "58.52", "2.91", "0.56"]
        rows = Dict(row["name"] => row for row in summary["rows"])
        @test (format_kg(rows["Water"]["kg"]), format_mol(rows["Water"]["mol"])) == ("8.500", "471.8")
        @test format_kg(rows["Ash"]["kg"]) == "0.127"
        @test format_kg(rows["Ash-free solids"]["kg"]) == "1.373"
        @test (format_kg(rows["C"]["kg"]), format_mol(rows["C"]["mol"])) == ("0.456", "37.96")
        @test isapprox(rows["Water"]["kg"], 8.5; rtol=1e-9) && isapprox(rows["Water"]["mol"], 471.8; rtol=1e-3)
        @test entered_sum_notice(summary["entered_sum"]) ==
              "The element percentages add up to 91.56 %. They have been scaled to 100 % (the Used column)."
        @test entered_sum_notice(100.0) === nothing && entered_sum_notice(99.999) === nothing
        @test format_mol(0.0) == "0" && format_mol(1234.4) == "1234"
        @test format_mol(999.96) == "1000" && format_mol(0.1591) == "0.1591" && format_mol(25.31) == "25.31"
        @test format_mol(1.0) == "1.000" && format_mol(0.0005) isa String
    end

    @testset "feed problems name their field" begin
        feed = merge(FEED, OrderedDict{String,Any}("ash" => OrderedDict{String,Any}("value" => 2, "unit" => "kg")))
        @test check_feed(feed)[2:3] == ("ash", "Ash can't be more than the solids.")
        @test check_feed(merge(FEED, OrderedDict{String,Any}("total_kg" => 0)))[2] == "total"
        @test check_feed(merge(FEED, OrderedDict{String,Any}("elements_wt_percent" => OrderedDict{String,Any}("C" => -1))))[2] == "element:C"
        @test check_feed(merge(FEED, OrderedDict{String,Any}("elements_wt_percent" => OrderedDict{String,Any}())))[2] == "elements"
    end

    @testset "conflict notice: the three conditions" begin
        data() = OrderedDict{String,Any}(
            "conditions" => OrderedDict{String,Any}("temperature_C" => Any[300, 400], "pressure_bar" => 240),
            "aqueous" => OrderedDict{String,Any}("cutoff_C" => 390.0),
            "constraints" => OrderedDict{String,Any}("carbon" => OrderedDict{String,Any}("aqueous" => OrderedDict{String,Any}("exact" => 0.03))))
        @test conflict_notice(data()) == "Aqueous carbon is set to exactly 3 % (Carbon split). At 390 °C and above there is no " *
                                         "aqueous phase, so that setting is ignored there. Below 390 °C it applies."
        @test occursin("At 663.15 K and above", conflict_notice(data(), "K"))
        d = data(); d["constraints"]["carbon"]["aqueous"] = OrderedDict{String,Any}("min" => 0.05)
        @test startswith(conflict_notice(d), "Aqueous carbon is set to at least 5 %")
        for change in (d -> d["aqueous"] = OrderedDict{String,Any}(),                                          # no cut-off
                       d -> d["constraints"]["carbon"]["aqueous"] = OrderedDict{String,Any}("max" => 0.05),    # not exactly / at least
                       d -> d["constraints"]["carbon"]["aqueous"] = OrderedDict{String,Any}("exact" => 0.0),   # zero
                       d -> d["conditions"]["temperature_C"] = Any[300, 380])                                  # all below
            d = data(); change(d)
            @test conflict_notice(d) === nothing
        end
    end

    @testset "initial guesses tooltip: singular and plural" begin
        @test startswith(initial_guesses_tooltip(1), "You will start from 1 initial guess. The optimization")
        @test startswith(initial_guesses_tooltip(5), "You will start from 5 initial guesses. The optimization")
        @test endswith(initial_guesses_tooltip(5), "they are listed on the Results tab so you can compare them.")
    end

    @testset "small helpers" begin
        c(pairs...) = OrderedDict{String,Any}(pairs...)
        @test point_count(c("temperature_C" => [300, 350, 400], "pressure_bar" => [240, 300])) == 6
        @test point_count(c("temperature_C" => [300, 400], "pressure_bar" => 240), 390) == 3
        @test point_count(c("points" => [[300, 240], [350, 250]])) == 2
        @test point_count(c("temperature_C" => [300])) === nothing
        @test (point_count_text(13), point_count_text(1), point_count_text(nothing)) == ("13 points", "1 point", "")
        @test element_error(element_symbol("cl")) === nothing
        @test element_error(element_symbol("qq")) == "\"Qq\" is not an element symbol."
        @test element_error("C", ["C"]) == "C is already in the table."
        @test parse_number(" 2.5 ") == 2.5 && parse_number("") === nothing
        @test_throws FieldError parse_number(""; required=true)
        @test (try parse_number("1.5", "Initial guesses"; integer=true) catch e; e.msg end) == "Initial guesses: give a whole number."
        @test problem_text("3 · Ash", "Ash can't be more than the solids.") == "Ash can't be more than the solids."
        @test problem_text("1 · Total mass [kg]", "Give a value.") == "Total mass: give a value."
        @test parse_points("573.15, 24\n623.15 30", "K", "MPa") == [[300, 240], [350, 300]]
        @test file_status(nothing, true) == "New case · not saved yet"
        @test file_status(joinpath(REPO, "cases", "x.toml"), true, REPO) == "cases/x.toml · unsaved changes"
        @test file_status(joinpath(REPO, "cases", "x.toml"), false, REPO) == "cases/x.toml"
        @test with_extension("a/b") == "a/b.toml" && with_extension("a/b.txt") == "a/b.txt" && with_extension("") == ""
    end

    @testset "SUPCRT92 band notice" begin
        d = HTGApp.new_document().data
        d["conditions"] = OrderedDict{String,Any}("temperature_C" => Any[300.0, 379.0, 386.0, 395.0], "pressure_bar" => 240.0)
        d["aqueous"] = OrderedDict{String,Any}("cutoff_C" => 390.0)
        # only electrolyte + SUPCRT92 + an aqueous phase; not at or above the cut-off
        @test band_points(d) == [(300.0, 240.0), (379.0, 240.0), (386.0, 240.0)]
        d["aqueous"]["water_model"] = "IAPWS-95"
        @test isempty(band_points(d))
        d["aqueous"]["water_model"] = "SUPCRT92"
        d["aqueous"]["model"] = "molecular"
        @test isempty(band_points(d))
        delete!(d["aqueous"], "model")
        d["aqueous"]["enabled"] = false
        @test isempty(band_points(d))
        delete!(d["aqueous"], "enabled")
        # the real check (SUPCRT92 at 240 bar: 379 °C in a band, 386 °C within 0.5 °C of one)
        hits = band_hits(band_points(d))
        @test hits == [(379.0, 240.0), (386.0, 240.0)]
        @test band_notice(hits) == "SUPCRT92 can fail at 379, 386 °C at 240 bar: in or within 0.5 °C of a band where " *
                                   "its water density does not converge. Results at these points may be wrong. Move " *
                                   "these temperatures, or choose IAPWS-95 as the water model."
        @test startswith(band_notice(hits, "K", "MPa"), "SUPCRT92 can fail at 652.15, 659.15 K at 24 MPa: in or within 0.5 K")
        @test band_notice([]) === nothing
        many = [(300.0 + i, P) for P in (240.0, 250.0) for i in 1:6]
        @test occursin("; 301, 302 °C at 250 bar and 4 more:", band_notice(many))
    end

    @testset "documents (Python: CaseDocument)" begin
        dir = mktempdir()
        # every example case opens and saves back to an equal case
        for name in readdir(joinpath(REPO, "cases"))
            endswith(name, ".toml") || continue
            doc = open_document(joinpath(REPO, "cases", name))
            @test doc.legacy_blocker === nothing && !dirty(doc)
            path = joinpath(dir, name)
            save_document(doc, path)
            @test HTGEquilibrium.dump_case(load_case(path)) == HTGEquilibrium.dump_case(load_case(joinpath(REPO, "cases", name)))
            again = open_document(path)   # a fresh file: the defaults are written out
            @test !dirty(again) && HTGEquilibrium.dump_case(document_case(again)) == HTGEquilibrium.dump_case(load_case(path))
        end
        # salts are kept, and noticed
        doc = open_document(joinpath(REPO, "cases", "manure_ions_minerals.toml"))
        @test has_salts(doc) && SALTS_NOTICE in notices(doc)
        path = save_document(doc, joinpath(dir, "salts.toml"))
        @test occursin("[feed.salts.", read(path, String))
        @test keys(load_case(path).feed.salts) == keys(load_case(joinpath(REPO, "cases", "manure_ions_minerals.toml")).feed.salts)
        # an old file is converted when opened; saving writes format 2
        old = joinpath(dir, "old.toml")
        cp(joinpath(REPO, "tests", "data", "legacy_cases", "manure_gas_min.toml"), old)
        doc = open_document(old)
        @test doc.converted_from_legacy && notices(doc)[1] == LEGACY_NOTICE
        save_document(doc)
        @test occursin("format = 2", read(old, String)) && !doc.converted_from_legacy && notices(doc) == [SALTS_NOTICE]
        # a dropped solver (owner decision PR #10): IPOPT, with the core's note
        text = replace(read(joinpath(REPO, "cases", "manure_gas_min.toml"), String), "method = \"ipopt\"" => "method = \"slsqp-fd\"\nfd_eps = 1e-8")
        doc = CaseDocument(text)
        @test !haskey(doc.data["solver"], "method") && !haskey(doc.data["solver"], "fd_eps")
        @test notices(doc)[1:2] == ["The solver 'slsqp-fd' is no longer available; this case now uses 'ipopt' (the default).",
                               "[solver] fd_eps was a setting of the removed solver 'slsqp-fd'; ignored."]
        @test document_case(doc).solver == "ipopt"
        # a file the form cannot show: kept as it is
        blocked = joinpath(dir, "unmodelled.toml")
        write(blocked, "name = \"old\"\n[conditions]\ntemperature_C = 300\npressure_bar = 240\n" *
                       "[feed]\nelements_wt_percent = { C = 50, H = 6, O = 40 }\nunmodelled_wt_percent = 4\n" *
                       "water_g_per_g_dry = 4\n[phases]\ngas = [\"water\", \"methane\", \"carbon dioxide\"]\n")
        doc = open_document(blocked)
        @test doc.legacy_blocker == "This case uses unmodelled_wt_percent, which the form can't show. Edit the file directly."
        @test_throws FieldError save_document(new_document())
    end

    config = mktempdir()
    withenv("HTG_CONFIG_DIR" => config) do
        state = State(; results_dir=mktempdir(), cases_dir=joinpath(REPO, "cases"))

        @testset "form: feed card numbers, ash error, elements" begin
            ed = editor(state)
            for (key, text) in (("feed.total", "10"), ("feed.solids", "15"), ("feed.ash", "8.45"), ("feed.element:C", "30.4"),
                                ("feed.element:H", "4.41"), ("feed.element:O", "53.58"), ("feed.element:N", "2.66"),
                                ("feed.element:S", "0.51"))
                edit!(ed, key, text)
            end
            v = form_view(ed)
            @test [v["feed"]["used"][e] for e in ("C", "H", "O", "N", "S")] == ["33.20", "4.82", "58.52", "2.91", "0.56"]
            @test v["feed"]["entered_sum"] == "91.56"
            @test v["feed"]["scale_notice"] == "The element percentages add up to 91.56 %. They have been scaled to 100 % (the Used column)."
            rows = Dict(r[1] => r[2:3] for r in v["feed"]["summary"])
            @test rows["Water"] == ["8.500", "471.8"] && rows["Ash"][1] == "0.127" && rows["Ash-free solids"][1] == "1.373"
            @test rows["C"] == ["0.456", "37.96"]
            feed = state.document.data["feed"]
            @test feed["total_kg"] == 10 && feed["solids"] == Dict("value" => 15, "unit" => "%")
            @test feed["ash"] == Dict("value" => 8.45, "unit" => "%") && v["ok"]
            # ash as kg, more than the 1.5 kg of solids
            share_unit!(ed, "ash", "kg")
            @test feed["ash"]["unit"] == "kg" && feed["ash"]["value"] ≈ 0.12675
            @test texts(ed)["feed.ash"] == "0.12675" && ed.texts_changed
            edit!(ed, "feed.ash", "2")
            v = form_view(ed)
            @test v["errors"]["feed.ash"] == "Ash can't be more than the solids."
            @test !v["ok"] && v["status"] == "✗ Ash can't be more than the solids."
            @test queue!(ed) === nothing && last(ed.notifications) == ("Not queued: Ash can't be more than the solids.", "negative")
            edit!(ed, "feed.ash", "0.1")
            @test !haskey(errors(ed), "feed.ash") && form_view(ed)["ok"]
            # adding and removing elements
            add_element!(ed, "qq")
            @test form_view(ed)["new_element_error"] == "\"Qq\" is not an element symbol."
            add_element!(ed, "cl")
            @test collect(keys(feed["elements_wt_percent"])) == ["C", "H", "O", "N", "S", "Cl"]
            @test form_view(ed)["new_element_error"] === nothing
            remove_element!(ed, "S")
            @test !haskey(feed["elements_wt_percent"], "S") && form_view(ed)["feed"]["elements"] == ["C", "H", "O", "N", "Cl"]
            # a field that is not a number blocks saving and names the field
            edit!(ed, "feed.total", "ten")
            v = form_view(ed)
            @test v["status"] == "✗ Total mass: \"ten\" is not a number." && feed["total_kg"] == 10 && !v["ok"]
            @test !write!(ed, joinpath(mktempdir(), "x.toml"))
        end

        @testset "form: conflict notice" begin
            ed = editor(state)
            data = state.document.data
            conflict() = form_view(ed)["aqueous"]["conflict"]
            edit!(ed, "conditions.temperatures", "300, 400")
            edit!(ed, "aqueous.cutoff", "390")
            carbon_choice!(ed, "aqueous", "exact")
            @test conflict() === nothing
            edit!(ed, "carbon.aqueous", "3")
            @test data["constraints"]["carbon"]["aqueous"] == Dict("exact" => 0.03)
            @test startswith(conflict(), "Aqueous carbon is set to exactly 3 % (Carbon split). At 390 °C and above")
            edit!(ed, "aqueous.cutoff", "")        # 1: the cut-off
            @test conflict() === nothing
            edit!(ed, "aqueous.cutoff", "390")
            @test conflict() !== nothing
            carbon_choice!(ed, "aqueous", "max")   # 2: the aqueous carbon choice
            @test conflict() === nothing && data["constraints"]["carbon"]["aqueous"] == Dict("max" => 0.03)
            carbon_choice!(ed, "aqueous", "min")
            @test startswith(conflict(), "Aqueous carbon is set to at least 3 %")
            edit!(ed, "carbon.aqueous", "0")
            @test conflict() === nothing
            edit!(ed, "carbon.aqueous", "3")
            edit!(ed, "conditions.temperatures", "300, 380")   # 3: a temperature at or above the cut-off
            @test conflict() === nothing
            edit!(ed, "conditions.temperatures", "300:400:50")
            @test conflict() !== nothing
            # the notice and the fields follow the temperature unit (a preference)
            set_unit!(ed, "temperature_unit", "K")
            @test occursin("At 663.15 K and above", conflict())
            @test texts(ed)["conditions.temperatures"] == "573.15, 623.15, 673.15" && ed.texts_changed
            @test state.settings["temperature_unit"] == "K"
            @test form_view(ed)["labels"]["conditions.temperatures"] == "Temperatures [K]"
            @test form_view(ed)["labels"]["aqueous.cutoff"] == "Remove at and above [K]"
            edit!(ed, "conditions.temperatures", "573.15, 773.15")
            @test data["conditions"]["temperature_C"] == [300, 500]
            set_unit!(ed, "temperature_unit", "C")
        end

        @testset "form: initial guesses" begin
            ed = editor(state)
            edit!(ed, "solver.starts", "1")
            @test startswith(form_view(ed)["solver"]["starts_tooltip"], "You will start from 1 initial guess.")
            edit!(ed, "solver.starts", "7")
            @test startswith(form_view(ed)["solver"]["starts_tooltip"], "You will start from 7 initial guesses.")
            @test state.document.data["solver"]["starts"] == 7
            edit!(ed, "solver.starts", "0")
            @test form_view(ed)["status"] == "✗ Initial guesses: give a value of 1 or more."
            @test state.document.data["solver"]["starts"] == 7
        end

        @testset "form: conditions and point count" begin
            ed = editor(state)
            data = state.document.data
            count() = form_view(ed)["conditions"]["count"]
            edit!(ed, "conditions.temperatures", "300:600 n=13")
            @test count() == "14 points"   # the template's cut-off (390 °C) adds a point
            edit!(ed, "aqueous.cutoff", "")
            @test count() == "13 points"
            set_unit!(ed, "pressure_unit", "MPa")
            edit!(ed, "conditions.pressures", "24, 30")
            @test data["conditions"] == Dict("temperature_C" => [300 + 25 * i for i in 0:12], "pressure_bar" => [240, 300])
            @test count() == "26 points"
            switch_conditions!(ed, "points")
            @test length(data["conditions"]["points"]) == 26 && !haskey(data["conditions"], "temperature_C")
            @test form_view(ed)["conditions"]["mode"] == "points"
            edit!(ed, "conditions.points", "300, 24\n350, 25")
            @test data["conditions"] == Dict("points" => [[300, 240], [350, 250]]) && count() == "2 points"
            set_input_formats!(ed, true)
            @test state.settings["input_formats_open"] === true && form_view(ed)["conditions"]["formats_open"]
            settings = HTGEquilibrium.load_settings(config)
            @test settings["pressure_unit"] == "MPa" && settings["input_formats_open"] === true
            set_unit!(ed, "pressure_unit", "bar")
            switch_conditions!(ed, "grid")
            @test data["conditions"] == Dict("temperature_C" => [300, 350], "pressure_bar" => [240, 250])
        end

        @testset "form: SUPCRT92 band notice (warn only)" begin
            calls = Ref(0)
            ed = editor(state; band_check=fake_bands([379.0]; calls))
            @test form_view(ed)["conditions"]["band"] === nothing
            n = calls[]
            edit!(ed, "conditions.temperatures", "300, 379")
            @test startswith(form_view(ed)["conditions"]["band"], "SUPCRT92 can fail at 379 °C at 240 bar")
            @test form_view(ed)["ok"]   # a warning: the case can still be saved and queued
            edit!(ed, "solver.starts", "2")   # the points did not change: no new check
            @test calls[] == n + 1
            set_water_model!(ed, "IAPWS-95")
            @test state.document.data["aqueous"]["water_model"] == "IAPWS-95" && form_view(ed)["conditions"]["band"] === nothing
            set_water_model!(ed, "SUPCRT92")
            @test !haskey(state.document.data["aqueous"], "water_model") && form_view(ed)["conditions"]["band"] !== nothing
            # the water model can be chosen only with the Electrolyte model and an aqueous phase
            @test form_view(ed)["aqueous"]["water_model_reason"] === nothing
            set_choice!(ed, "aqueous_model", "molecular")
            @test form_view(ed)["aqueous"]["water_model_reason"] == WATER_MODEL_MOLECULAR
            @test form_view(ed)["conditions"]["band"] === nothing
            set_choice!(ed, "aqueous_model", "electrolyte")
            set_choice!(ed, "aqueous_enabled", false)
            @test state.document.data["aqueous"]["enabled"] === false
            @test form_view(ed)["aqueous"]["water_model_reason"] == WATER_MODEL_NO_AQUEOUS
            set_choice!(ed, "aqueous_enabled", true)
            @test !haskey(state.document.data["aqueous"], "enabled")
            # in the background: the result comes later, through on_change
            changed = Channel{Bool}(1)
            ed2 = CaseEditor(state; band_check=fake_bands([386.0]), on_change=() -> put!(changed, true))
            edit!(ed2, "conditions.temperatures", "386")
            @test take!(changed) && startswith(form_view(ed2)["conditions"]["band"], "SUPCRT92 can fail at 386 °C")
        end

        @testset "form: simple choices write defaults as no key" begin
            ed = editor(state)
            d = state.document.data
            set_choice!(ed, "method", "slsqp")
            @test d["solver"]["method"] == "slsqp" && form_view(ed)["solver"]["method"] == "slsqp"
            set_choice!(ed, "continuation", false)
            @test d["sweep"]["continuation"] === false && !form_view(ed)["solver"]["continuation"]
            set_choice!(ed, "continuation", true)
            @test HTGApp.getpath(d, ("sweep", "continuation")) === nothing
            set_choice!(ed, "stability", "off")
            @test d["solver"]["stability_test"] === false && form_view(ed)["advanced"]["stability"] == "off"
            set_choice!(ed, "stability", "default")
            @test !haskey(d["solver"], "stability_test")
            set_choice!(ed, "solid_target", "all_solids")
            @test d["constraints"]["solid_carbon_target"] == "all_solids"
            set_choice!(ed, "solid_target", "char")
            @test HTGApp.getpath(d, ("constraints", "solid_carbon_target")) === nothing
            set_choice!(ed, "initialization", "lp")
            @test d["initial_guess"]["initialization"] == "lp"
            @test document_case(state.document).solver == "slsqp"
        end

        model = load_model(; user_data=nothing)

        @testset "form: species amounts" begin
            ed = editor(state)
            switch_feed!(ed, "species")
            @test form_view(ed)["feed"]["summary_message"] == LOADING_MESSAGE && !form_view(ed)["ok"]
            use_model!(state, model)
            HTGApp.refresh!(ed)
            feed = state.document.data["feed"]
            @test feed == Dict("species" => Dict("water" => Dict("value" => 4.5, "unit" => "kg")))   # no feed species before
            switch_feed!(ed, "elements"); delete!(ed.feed_stash, "species"); switch_feed!(ed, "species")
            feed = state.document.data["feed"]
            @test feed == Dict("species" => Dict("water" => Dict("value" => 4.5, "unit" => "kg"),
                                                 "cellulose" => Dict("value" => 0.5, "unit" => "kg")))
            @test form_view(ed)["feed"]["caption"] == HTGApp.FEED_CAPTION_SPECIES
            edit!(ed, "feed.water", "9")
            @test feed["species"]["water"] == Dict("value" => 9, "unit" => "kg")
            water_unit!(ed, "mol")
            @test feed["species"]["water"]["unit"] == "mol"
            @test feed["species"]["water"]["value"] ≈ 9000 / 18.015 rtol = 1e-5
            add_species!(ed)
            @test form_view(ed)["status"] == "✗ Species: choose a species or remove the row."
            species_name!(ed, 1, "cellulose")
            @test errors(ed)["feed.species-name:1"] == "cellulose is already in the feed."
            species_name!(ed, 1, "methanol")
            @test collect(keys(state.document.data["feed"]["species"])) == ["water", "cellulose", "methanol"]
            @test form_view(ed)["ok"]
            @test case_from_data(state.document.data).feed_v2["species"]["methanol"] == Dict("value" => 1.0, "unit" => "kg")
            rows = Dict(r[1] => r for r in form_view(ed)["feed"]["summary"])
            @test haskey(rows, "Total")
            switch_feed!(ed, "elements")
            @test haskey(state.document.data["feed"], "total_kg")
            switch_feed!(ed, "species")   # the species come back
            @test haskey(state.document.data["feed"]["species"], "methanol")
            remove_species!(ed, 1)
            @test collect(keys(state.document.data["feed"]["species"])) == ["water", "cellulose"]
        end

        @testset "form: open, save writes format 2, reload keeps unsaved changes" begin
            dir = mktempdir()
            path = joinpath(dir, "old.toml")
            cp(joinpath(REPO, "tests", "data", "legacy_cases", "manure_gas_min.toml"), path)
            ed = editor(state)
            @test open!(ed, path)
            v = form_view(ed)
            @test LEGACY_NOTICE in v["notices"] && ed.texts_changed
            edit!(ed, "solver.starts", "3")
            @test endswith(form_view(ed)["file"]["status"], "old.toml · unsaved changes")
            # a reloaded page: a new editor over the same document still knows
            ed2 = CaseEditor(state; async_bands=false, band_check=fake_bands(Float64[]))
            @test ed2.texts["solver.starts"] == "3" && endswith(form_view(ed2)["file"]["status"], "· unsaved changes")
            @test save!(ed)
            @test !dirty(state.document) && !(LEGACY_NOTICE in form_view(ed)["notices"])
            @test last(ed.notifications) == ("Saved $(realpath(path))", "positive")
            text = read(path, String)
            @test occursin("format = 2", text) && occursin("starts = 3", text) && !occursin("multistart", text)
            @test load_case(path).starts == 3
            # Save as
            other = joinpath(dir, "sub", "copy.toml")
            @test write!(ed, other) && state.document.path == realpath(other)
            # a new case has no file: Save asks for one
            new!(ed)
            @test !save!(ed) && form_view(ed)["file"]["status"] == "New case · not saved yet"
            @test HTGApp.default_save_path(ed) == joinpath(state.cases_dir, "new_case.toml")
            @test !open!(ed, joinpath(dir, "missing.toml"))
            @test startswith(last(ed.notifications)[1], "Could not open ")
        end

        @testset "form: salts and blocked files" begin
            ed = editor(state)
            open!(ed, joinpath(REPO, "cases", "manure_ions_minerals.toml"))
            v = form_view(ed)
            @test v["salts"] && v["ok"] && !(SALTS_NOTICE in v["notices"])   # the page shows it above the Feed card
            blocked = joinpath(mktempdir(), "unmodelled.toml")
            write(blocked, "name = \"old\"\n[conditions]\ntemperature_C = 300\npressure_bar = 240\n" *
                           "[feed]\nelements_wt_percent = { C = 50, H = 6, O = 40 }\nunmodelled_wt_percent = 4\n" *
                           "water_g_per_g_dry = 4\n[phases]\ngas = [\"water\", \"methane\", \"carbon dioxide\"]\n")
            open!(ed, blocked)
            v = form_view(ed)
            @test v["blocked"] && isempty(ed.fields) && v["ok"]
            @test "This case uses unmodelled_wt_percent, which the form can't show. Edit the file directly." in v["notices"]
            @test v["status"] == "Ready to run: 1 point in 1 independent piece."
        end

        @testset "form: check" begin
            ed = editor(state)
            c = check_case(ed)
            check_done!(ed, c, preflight_messages(state.model, c))
            @test ed.check["open"] && ed.check["ok"] && ed.check["title"] == "All 5 points passed the checks."
            carbon_choice!(ed, "gas", "min"); edit!(ed, "carbon.gas", "90")
            carbon_choice!(ed, "solid", "min"); edit!(ed, "carbon.solid", "20")
            c = check_case(ed)
            messages = preflight_messages(state.model, c)
            check_done!(ed, c, messages)
            @test !ed.check["ok"] && ed.check["title"] == "Problems found"
            @test any(m -> occursin("The carbon shares add up to more than 100 %.", m), ed.check["messages"])
            loading = State(; settings=state.settings, cases_dir=state.cases_dir)
            ed = editor(loading)
            @test check_case(ed) === nothing && last(ed.notifications) == (LOADING_MESSAGE, "warning")
        end

        @testset "form: add to queue" begin
            ed = editor(state)
            @test queue!(ed) === nothing   # no pool yet: the species data are loading
            @test last(ed.notifications) == ("Not queued: " * HTGApp.LOADING_MESSAGE, "negative")
            pool = WorkerPool(1; model, warm_up=false, chain_solver=stand_in)
            state.pool = pool
            HTGEquilibrium.start(pool)
            try
                run = queue!(ed)
                @test run !== nothing && last(ed.notifications) == ("Queued form_test: 5 points", "positive")
                @test wait(run, 60)
            finally
                shutdown(pool)
                state.pool = nothing
            end
        end
    end

    @testset "page markup" begin
        html = case_tab_html()
        for mark in ("case-file", "case-new", "case-save", "case-save-as", "case-check", "case-queue", "case-status",
                     "case-conditions-mode", "case-temperature-unit", "case-point-count", "case-band", "case-input-formats",
                     "case-feed-mode", "case-feed-summary", "case-entered-sum", "case-add-element", "case-scale-notice",
                     "case-phase-gas", "case-aqueous-enabled", "case-water-model", "case-conflict", "case-method",
                     "case-starts-tooltip", "case-advanced", "case-check-dialog", "case-saveas-dialog")
            @test occursin("data-mark=\"$mark\"", html)
        end
        @test occursin(HTGApp._attr(HTGApp.WATER_MODEL_TOOLTIP), html) && occursin(HTGApp._attr(HTGApp.METHOD_TOOLTIP), html)
        @test !occursin("slsqp-fd", html) && !occursin("trust-constr", html)
    end
end
