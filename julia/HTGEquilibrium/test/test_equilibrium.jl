# Complete equilibrium problems: tests/reference/data/solves.json (Python:
# tests/test_reference_solves.py), solving checks of tests/test_backend.py, the
# Result data, and solving on threads.
#
# Two levels, as in Python:
# - setup: the same balances, feasible start and Gibbs energy function as the
#   original (G evaluated at the original's solution must agree to ~1e-10);
# - solve: the optimizer reaches the original's solution.
#
# Most reference points were solved with the original finite-difference SLSQP
# ("slsqp-fd"), which was dropped (owner, PR #10: "the stored reference used a bad
# optimizer", so these points are a rough guide only): they are solved with IPOPT
# (the default), as the Python test already does for fig78_380C. Their G must not
# be higher than the reference; where IPOPT finds a lower G (LOWER_THAN_FD) the point is checked
# against the Python app's own IPOPT solution (test/data/solves_python_ipopt.json,
# written by test/dump_solves.py) and must be stationary. Every point is also
# compared with that file.
#
# HTG_FULL_TESTS=1 solves every reference point; the default is a representative
# subset (DEFAULT_SOLVES) that covers each kind of check.
using Test
using HTGEquilibrium
using HTGEquilibrium: OrderedDict
using JSON
using LinearAlgebra: BLAS
isdefined(Main, :DataTestSetup) || include("data_setup.jl")
isdefined(Main, :SolvesTestSetup) || include("solves_setup.jl")
using Main.DataTestSetup: MODELS, REPO, reldev
using Main.SolvesTestSetup

BLAS.set_num_threads(1)   # M1: Ipopt/MUMPS run faster and threads do not oversubscribe

const SOLVE_TEMPLATE = joinpath(REPO, "src", "htg", "templates", "case_template.toml")
const PYTHON_IPOPT_FILE = joinpath(@__DIR__, "data", "solves_python_ipopt.json")
const PYTHON_IPOPT = isfile(PYTHON_IPOPT_FILE) ? JSON.parsefile(PYTHON_IPOPT_FILE) : nothing

# Python: IMPROVED. The package finds a lower G than the stored reference, and the
# reference case whose solution it must reach instead (the original IPOPT setting
# mu_init = 1e-6 stopped early at full_text 250 C).
const IMPROVED = Dict("full_text_250C_ipopt" => "full_text_250C")
# Python: FD_STOPS_EARLY. The original finite-difference SLSQP stops early at
# fig78 380 C (in Python too); the default exact-gradient solver must reach the reference.
const FD_STOPS_EARLY = Set(["fig78_380C"])
# slsqp-fd reference points where IPOPT (Python's and Julia's alike) reaches a G
# lower than the reference by more than 1e-9 relative: the finite-difference
# SLSQP stopped slightly early. Checked against Python's IPOPT solution instead.
const LOWER_THAN_FD = Set(["cellulose_400C", "full_text_250C", "full_text_gas_min_300C", "configurable_300C"])
# Default subset: IPOPT, multistart (random starts), SLSQP, slsqp-fd points (molecular
# and electrolyte, carbon bounds, a species feed, a lower G), IMPROVED, FD_STOPS_EARLY.
const DEFAULT_SOLVES = Set(["fig78_380C", "fig78_450C_ipopt", "fig78_450C_ipopt_multistart",
                            "fig78_450C_slsqp_exact", "fig78_450C_allocations", "fig78_450C_all_solids",
                            "cellulose_species_feed_450C", "full_text_250C", "full_text_250C_ipopt",
                            "configurable_300C_ipopt", "configurable_450C"])

uses_fd(entry) = get(get(entry["result"], "optimizer", Dict()), "method", "slsqp-fd") == "slsqp-fd"
"""The reference case (slsqp-fd points with the default solver, IPOPT: see case_from_legacy)."""
solve_case(entry) = case_from_legacy(entry)
amount_values(triples) = collect(values(legacy_pairs(triples)))
"""numpy assert_allclose: |a - b| <= atol + rtol |b| elementwise; also the worst ratio."""
allclose(a, b; rtol, atol) = all(abs.(a .- b) .<= atol .+ rtol .* abs.(b))
worst(a, b; rtol, atol) = maximum(abs.(a .- b) ./ (atol .+ rtol .* abs.(b)))

const SUMMARY = OrderedDict{String,Any}()   # per point: deviations and iterations, printed at the end

@testset "equilibrium" begin
    @testset "setup and Gibbs function: $(entry["name"])" for entry in SOLVED
        ref = entry["result"]
        s = setup(MODELS[entry["model"]], case_from_legacy(entry))
        expected = legacy_pairs(ref["amounts"])
        @test s.pairs == collect(keys(expected))
        @test allclose(s.initial_amounts, amount_values(ref["initial_amounts"]); rtol=1e-9, atol=1e-11)
        @test s.elements == ref["elements"]
        @test allclose(s.b, Float64.(ref["feed_atoms"]); rtol=1e-12, atol=1e-14)
        @test reldev(s.C_mol_per_g, ref["feed_C_mol_per_g_dry"]) <= 1e-12
        terms = legacy_terms(s.phase_G(collect(values(expected))))
        for (phase, value) in ref["phase_gibbs_J_per_mol_feed_C"]
            @test isapprox(terms[phase], value; rtol=1e-10, atol=1e-6)
        end
        # the start is feasible and the LP part matches Python's (HiGHS 1.12.0, as SciPy 1.18)
        @test maximum(abs.(s.A_full * s.initial_amounts .- s.b_full)) <= 1e-8
        @test s.x_start == s.initial_amounts[s.free]
    end

    @testset "solve: $(entry["name"])" for entry in SOLVED
        name = entry["name"]
        (FULL_TESTS || name in DEFAULT_SOLVES) || continue
        ref = entry["result"]
        expected = legacy_pairs(ref["amounts"])
        seconds = @elapsed r = solve(MODELS[entry["model"]], solve_case(entry))
        x = collect(values(r.amounts))
        G = r.gibbs_J_per_mol_C
        info = SUMMARY[name] = OrderedDict{String,Any}("reference" => uses_fd(entry) ? "slsqp-fd" : r.optimizer["method"],
                                                         "dG_rel" => (G - ref["gibbs_J_per_mol_feed_C"]) / abs(ref["gibbs_J_per_mol_feed_C"]),
                                                         "seconds" => seconds)
        @test collect(keys(r.amounts)) == collect(keys(expected))
        @test r.balance_residual < 1e-7
        rules = ref["phase_allocation"]
        @test r.phase_rules["aqueous_disabled"] == rules["aqueous_disabled"]
        @test r.phase_rules["removed_aqueous_species"] == rules["removed_aqueous_species"]
        @test r.phase_rules["initial_composition_reset"] == rules["initial_composition_reset"]
        python = PYTHON_IPOPT === nothing ? nothing : PYTHON_IPOPT[name]
        if haskey(IMPROVED, name)
            better = REFERENCE[IMPROVED[name]]["result"]
            @test G < ref["gibbs_J_per_mol_feed_C"]
            @test reldev(G, better["gibbs_J_per_mol_feed_C"]) <= 1e-8
            # G is flat in the gas/aqueous water split here, so amounts agree less tightly than G.
            @test allclose(x, amount_values(better["amounts"]); rtol=1e-2, atol=1e-4)
            info["composition"] = worst(x, amount_values(better["amounts"]); rtol=1e-3, atol=1e-4)
        elseif uses_fd(entry)
            # IPOPT instead of the finite-difference SLSQP: never a higher G
            @test G <= ref["gibbs_J_per_mol_feed_C"] + 1e-9 * abs(ref["gibbs_J_per_mol_feed_C"])
            if name in LOWER_THAN_FD
                @test reldev(G, ref["gibbs_J_per_mol_feed_C"]) > 1e-9
                @test r.optimizer["stationarity_major_RT"] < POLISH_THRESHOLD_RT
                if python !== nothing   # the Python app's IPOPT reaches the same lower G
                    @test reldev(G, python["gibbs_J_per_mol_C"]) <= 1e-9
                    @test allclose(x, Float64[v for (_, _, v) in python["amounts"]]; rtol=1e-3, atol=1e-4)
                end
                info["composition_vs_reference"] = worst(x, collect(values(expected)); rtol=1e-3, atol=1e-4)
            else
                @test reldev(G, ref["gibbs_J_per_mol_feed_C"]) <= 1e-9
                # The aqueous phase is nearly empty at some points; trace amounts left in it barely change G.
                @test allclose(x, collect(values(expected)); rtol=1e-3, atol=1e-4)
                info["composition"] = worst(x, collect(values(expected)); rtol=1e-3, atol=1e-4)
            end
            @test r.optimizer["method"] == "ipopt"
        else
            r.optimizer["method"] == "ipopt" && @test !r.optimizer["polish"]["used"]
            @test reldev(G, ref["gibbs_J_per_mol_feed_C"]) <= 1e-9
            @test allclose(x, collect(values(expected)); rtol=1e-3, atol=1e-4)
            info["composition"] = worst(x, collect(values(expected)); rtol=1e-3, atol=1e-4)
            for (phase, value) in ref["phase_gibbs_J_per_mol_feed_C"]
                @test isapprox(legacy_terms(r.phase_gibbs_J_per_mol_C)[phase], value; rtol=1e-8, atol=1e-3)
            end
            @test isapprox(r.gas_carbon_fraction, ref["gas_carbon_fraction_actual"]; atol=1e-7)
            @test r.optimizer["method"] == ref["optimizer"]["method"]
            # Same starts reaching the same G (random starts: Julia's RNG gives other
            # points; they reach the same G here); iteration counts differ (exact Hessian).
            @test length(r.optimizer["starts"]) == length(ref["optimizer"]["starts"])
            for (new, old) in zip(r.optimizer["starts"], ref["optimizer"]["starts"])
                @test replace(new["start"], "_" => " ") == old["start"]   # phase labels renamed
                @test new["message"] == old["message"]
                @test reldev(new["G"], old["G"]) <= 1e-9
            end
        end
        if python !== nothing
            # the Python app's IPOPT (finite-difference Hessian) from the same starts
            @test reldev(G, python["gibbs_J_per_mol_C"]) <= 1e-9
            px = Float64[v for (_, _, v) in python["amounts"]]
            @test allclose(x, px; rtol=1e-3, atol=1e-4)
            info["dG_rel_python_ipopt"] = (G - python["gibbs_J_per_mol_C"]) / abs(python["gibbs_J_per_mol_C"])
            info["composition_vs_python_ipopt"] = worst(x, px; rtol=1e-3, atol=1e-4)
            if r.optimizer["method"] == python["solver"] == "ipopt"
                @test [s["message"] for s in r.optimizer["starts"]] == [s["message"] for s in python["starts"]]
                @test r.optimizer["polish"]["used"] == python["polish_used"]
            end
            info["nit"] = [get(s, "nit", nothing) for s in r.optimizer["starts"]]
            info["nit_python"] = [s["nit"] for s in python["starts"]]
        end
        # the result's own bookkeeping
        @test r.alternatives[1].label == "A" && r.alternatives[1].gibbs == G
        @test r.product_atoms ≈ r.feed_atoms rtol = 1e-9 atol = 1e-9
        @test abs(r.charge_residual) < 1e-7
        @test sum(values(r.phase_gibbs_J_per_mol_C)) ≈ G rtol = 1e-14
        @test r.water_model == "SUPCRT92"
        @test r.temperature_C == entry["settings"]["temperature_C"]
    end

    @testset "expected failure: $(entry["name"])" for entry in FAILED
        err = try
            solve(MODELS[entry["model"]], case_from_legacy(entry))
        catch e
            e
        end
        @test err isa SolveError
        # Same stage and cause; the advice at the end of some messages was reworded.
        expected = split(entry["error"], ": "; limit=2)[2]
        first_sentence = split(expected, ". ")[1]
        @test startswith(err.msg, first_sentence)
    end

    @testset "IPOPT polish safeguard" begin
        # With the original IPOPT setting (mu_init = 1e-6), Python's IPOPT stopped early at
        # full_text 250 C; the SLSQP polish then reached the lower solution.
        entry, better = REFERENCE["full_text_250C_ipopt"], REFERENCE["full_text_250C"]["result"]
        case = replace(case_from_legacy(entry); solver_options=Dict("ipopt_options" => Dict("mu_init" => 1e-6)))
        r = solve(MODELS[entry["model"]], case)
        polish = r.optimizer["polish"]
        SUMMARY["polish safeguard"] = OrderedDict{String,Any}(
            "used" => polish["used"], "stationarity_major_RT" => polish["stationarity_major_RT"],
            "G_ipopt_dev" => get(polish, "G_ipopt", NaN) - entry["result"]["gibbs_J_per_mol_feed_C"])
        @test polish["used"] && polish["stationarity_major_RT"] > 1
        @test reldev(polish["G_ipopt"], entry["result"]["gibbs_J_per_mol_feed_C"]) <= 1e-9
        @test reldev(r.gibbs_J_per_mol_C, better["gibbs_J_per_mol_feed_C"]) <= 1e-8
        # a polished result: its iterations add SLSQP's evaluations to IPOPT's
        @test r.iterations > r.optimizer["starts"][1]["nit"]
    end

    @testset "a case asking for slsqp-fd or trust-constr is solved with ipopt (owner, PR #10)" begin
        entry = REFERENCE["fig78_450C_ipopt"]
        model, case = MODELS[entry["model"]], case_from_legacy(entry)
        want = solve(model, case)
        for solver in ("slsqp-fd", "trust-constr")
            old = replace(case; solver)
            @test old.solver == "ipopt" && length(old.legacy_notes) == 1
            r = solve(model, old)
            @test r.optimizer["method"] == "ipopt" && isapprox(r.gibbs_J_per_mol_C, want.gibbs_J_per_mol_C; rtol=1e-12)
        end
    end

    @testset "points of a sweep case (Model.solve)" begin
        model = MODELS["base"]
        c = replace(load_case(SOLVE_TEMPLATE); temperature_C=[350, 400], pressure_bar=240,
                    carbon_bounds=OrderedDict("aqueous" => OrderedDict("exact" => 0.03)), aqueous_cutoff_C=390,
                    continuation=false)
        err = try solve(model, c) catch e e end
        @test err isa EquilibriumError && startswith(err.msg, "This case has several points;")
        @test_throws EquilibriumError setup(model, c)   # a single temperature and pressure is needed
        # aqueous = {exact = 0.03} with the cut-off inside the range: every point solves,
        # and the points at and above the cut-off say that the bound was dropped (test_backend.py).
        results = [solve(model, c, T) for T in temperatures(c)]
        @test [r.temperature_C for r in results] == [350.0, 390.0, 400.0]
        @test [notes(r) for r in results] == ["", "aqueous carbon bound dropped", "aqueous carbon bound dropped"]
        below = results[1]
        aqueous_carbon = sum(v * get(below.species_data[k]["elements"], "C", 0) for (k, v) in below.amounts
                             if k[1] == "aqueous")
        @test aqueous_carbon / below.feed_atoms[findfirst(==("C"), below.elements)] ≈ 0.03 atol = 1e-8
        @test all(r -> r.carbon_bound_violation < 1e-7, results)
        # several pressures need both
        c2 = replace(c; pressure_bar=[240, 300])
        err = try solve(model, c2, 400) catch e e end
        @test err isa CaseError && err.msg == "This case has several pressures; give pressure_bar too."
    end

    @testset "general carbon bounds are respected" begin
        c = replace(load_case(SOLVE_TEMPLATE); temperature_C=450,
                    carbon_bounds=OrderedDict("gas" => OrderedDict("min" => 0.5, "max" => 0.9),
                                              "solid" => OrderedDict("max" => 0.02)), aqueous_cutoff_C=390)
        r = solve(MODELS["base"], c)
        feed_C = r.feed_atoms[findfirst(==("C"), r.elements)]
        share = Dict{String,Float64}()
        for ((phase, s), v) in r.amounts
            share[phase] = get(share, phase, 0.0) + v * get(r.species_data[(phase, s)]["elements"], "C", 0) / feed_C
        end
        @test 0.5 - 1e-7 <= share["gas"] <= 0.9 + 1e-7
        @test share["solid"] <= 0.02 + 1e-7
        @test r.carbon_bound_violation < 1e-7
    end

    @testset "Result data" begin
        entry = REFERENCE["fig78_450C_ipopt"]
        r = solve(MODELS[entry["model"]], solve_case(entry))
        @test n_results(r) == 1 && notes(r) == ""
        @test for_alternative(r, "A").amounts == r.amounts
        @test_throws KeyError for_alternative(r, "B")
        # a second alternative: amounts, G, gas carbon and phase G follow it
        alt_amounts = copy(r.amounts)
        pair = first(k for k in keys(alt_amounts) if k[1] == "gas" && get(r.species_data[k]["elements"], "C", 0) > 0)
        alt_amounts[pair] *= 0.5
        push!(r.alternatives, Alternative("B", r.gibbs_J_per_mol_C + 10, alt_amounts, 1, ["random 0"];
                                          phase_gibbs=Dict("gas" => -1.0)))
        b = for_alternative(r, "B")
        @test n_results(r) == 2 && b.gibbs_J_per_mol_C == r.gibbs_J_per_mol_C + 10
        @test b.amounts[pair] == r.amounts[pair] / 2 && b.gas_carbon_fraction < r.gas_carbon_fraction
        @test b.phase_gibbs_J_per_mol_C["gas"] == -1.0 && b.phase_gibbs_J_per_mol_C["solid"] == 0.0
        @test r.amounts[pair] == alt_amounts[pair] * 2   # the original is unchanged
        @test all(haskey(d, "formula") && haskey(d, "elements") && haskey(d, "charge")
                                    for d in values(r.species_data))
        text = sprint(show, MIME("text/plain"), r)
        @test startswith(text, "Result at 450 °C, 240 bar: G = ") || startswith(text, "Result at 450 °C")
        @test occursin("mol/mol C", text)
    end

    @testset "threads: solves on $(Threads.nthreads()) threads equal the serial ones" begin
        names = ["fig78_450C_ipopt", "fig78_380C", "fig78_450C_allocations", "configurable_300C_ipopt",
                 "cellulose_species_feed_450C", "full_text_250C_ipopt"]
        FULL_TESTS && append!(names, ["fig78_450C_ipopt_multistart", "configurable_450C", "fig78_450C_slsqp_exact"])
        solve_named(name) = solve(MODELS[REFERENCE[name]["model"]], solve_case(REFERENCE[name]))
        serial = Dict(name => solve_named(name) for name in names)
        tasks = [(name, Threads.@spawn solve_named(name)) for _ in 1:4 for name in names]
        same = 0
        for (name, task) in tasks
            r, s = fetch(task), serial[name]
            ok = r.amounts == s.amounts && r.gibbs_J_per_mol_C == s.gibbs_J_per_mol_C &&
                 r.iterations == s.iterations && r.optimizer["nfev"] == s.optimizer["nfev"]
            same += ok
            @test ok
        end
        SUMMARY["threads"] = OrderedDict{String,Any}("threads" => Threads.nthreads(), "solves" => length(tasks),
                                                     "bitwise_equal" => same)
        Threads.nthreads() == 1 && @info "threads test ran on 1 thread; set JULIA_NUM_THREADS=4"
    end
end

let lines = ["solver summary ($(FULL_TESTS ? "HTG_FULL_TESTS" : "default subset")):"]
    for (name, info) in SUMMARY
        push!(lines, "  " * rpad(name, 30) * join(("$k=$(v isa AbstractFloat ? round(v; sigdigits=3) : v)"
                                                   for (k, v) in info), " "))
    end
    @info join(lines, "\n")
end
