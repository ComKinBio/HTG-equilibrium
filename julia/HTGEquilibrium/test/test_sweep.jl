# Sweeps and preflight checks (Python: htg/sweep.py, tests/test_grid.py, Model.check).
#
# Besides the Python tests, the solutions of the sweeps are checked on their
# own (owner, PR #10: tests that do not depend on the old Python code): element
# balances close, amounts are non-negative, the aqueous phase is electroneutral,
# carbon shares add up, and warm starts along a chain reach the same minimum as
# cold starts.
using Test
using HTGEquilibrium
using HTGEquilibrium: OrderedDict, save
using LinearAlgebra: BLAS
isdefined(Main, :DataTestSetup) || include("data_setup.jl")
using Main.DataTestSetup: MODELS, REPO

const SWEEP_TEMPLATE = joinpath(REPO, "src", "htg", "templates", "case_template.toml")
const FULL_SWEEP_TESTS = get(ENV, "HTG_FULL_TESTS", "") in ("1", "true", "yes")
template_case(; changes...) = replace(load_case(SWEEP_TEMPLATE); changes...)

"""Checks of a solution that need no reference: balances, signs, charge, carbon shares."""
function check_solution(r::Result)
    ok = true
    ok &= all(isapprox.(r.product_atoms, r.feed_atoms; rtol=1e-9, atol=1e-9))   # element balances
    ok &= r.balance_residual <= 1e-7
    ok &= all(>=(0), values(r.amounts))                                      # amounts non-negative
    ok &= abs(r.charge_residual) <= 1e-9                                     # electroneutral aqueous phase
    ok &= 0 <= r.gas_carbon_fraction <= 1 + 1e-9
    phases = phase_summary(r)
    feed_C = r.feed_atoms[findfirst(==("C"), r.elements)]
    product_C = r.product_atoms[findfirst(==("C"), r.elements)]
    ok &= isapprox(sum(phases["carbon_percent"]), 100 * product_C / feed_C; rtol=1e-9)
    ok &= isapprox(sum(phases["gibbs_J_per_mol_C"]), r.gibbs_J_per_mol_C; rtol=1e-9)
    species = species_table(r)
    ok &= isapprox(sum(species["mol_per_mol_C"]), sum(phases["mol_per_mol_C"]); rtol=1e-12)
    ok || @info "solution check failed" r
    return ok
end

same_results(a::SweepResult, b::SweepResult) =
    length(a.points) == length(b.points) &&
    all(p[1:2] == q[1:2] && typeof(p[3]) == typeof(q[3]) &&
        (p[3] isa String ? p[3] == q[3] :
         p[3].amounts == q[3].amounts && p[3].gibbs_J_per_mol_C == q[3].gibbs_J_per_mol_C &&
         p[3].iterations == q[3].iterations)
        for (p, q) in zip(a.points, b.points))

@testset "sweep" begin
    BLAS.set_num_threads(1)
    model = MODELS["base"]

    @testset "grids, explicit points, chains (test_grid.py)" begin
        c = template_case(temperature_C=[400, 300], pressure_bar=[250, 240], aqueous_cutoff_C=nothing)
        @test grid(c) == [(300.0, 240.0), (400.0, 240.0), (300.0, 250.0), (400.0, 250.0)]
        @test temperatures(c) == [300.0, 400.0] && pressures(c) == [240.0, 250.0]
        @test chains(c) == [[(300.0, 240.0), (400.0, 240.0)], [(300.0, 250.0), (400.0, 250.0)]]
        @test chains(replace(c; continuation=false)) == [[p] for p in grid(c)]
        e = template_case(temperature_C=nothing, pressure_bar=nothing, points=[[450, 300], [300, 240], [450, 300]])
        @test grid(e) == [(450.0, 300.0), (300.0, 240.0)] && chains(e) == [[(450.0, 300.0), (300.0, 240.0)]]
        r = template_case(temperature_C=[300, 400], pressure_bar=240, aqueous_cutoff_C=385,
                          refine=OrderedDict("start_C" => 350, "stop_C" => 370, "step_C" => 10))
        @test temperatures(r) == [300.0, 350.0, 360.0, 370.0, 385.0, 400.0]
        m = template_case(temperature_C=[300], pressure_bar=[240, 250])
        @test_throws CaseError at(m, 300)
        point = at(m, 300, 250)
        @test (point.temperature_C, point.pressure_bar, point.points) == (300.0, 250.0, nothing)
        @test_throws CaseError validate(template_case(temperature_C=nothing, pressure_bar=nothing, points=[[300]]))
        dir = mktempdir()
        for x in (template_case(pressure_bar=[240.0, 300.0]),
                  template_case(temperature_C=nothing, pressure_bar=nothing, points=[[300.0, 240.0], [400.0, 250.0]]))
            save(x, joinpath(dir, "c.toml"))
            @test load_case(joinpath(dir, "c.toml")) == x
        end
    end

    @testset "preflight (Model.check) with Python's messages" begin
        t = check(model, load_case(SWEEP_TEMPLATE))
        @test t.columns == ["T_C", "P_bar", "status", "aqueous_removed", "error"]
        @test t["T_C"] == [300.0, 350.0, 390.0, 400.0, 450.0] && all(==("ready"), t["status"])
        @test t["aqueous_removed"] == [false, false, true, true, true] && all(==(""), t["error"])
        c2 = template_case(aqueous_cutoff_C=nothing, temperature_C=[300, 400], pressure_bar=[240])
        message = "ValueError: Aqueous electrolyte model invalid here: water density 149 kg/m3 (needs >= 350). " *
                  "Set aqueous cutoff_C below this temperature."
        t2 = preflight(model, c2)
        @test t2["status"] == ["ready", "invalid"] && t2["error"] == ["", message]
        @test preflight_messages(model, c2) == ["400 °C, 240 bar: " * message]
        phases = OrderedDict{String,Any}(load_case(SWEEP_TEMPLATE).phases)
        c3 = template_case(temperature_C=[300], phases=merge(phases, Dict("gas" => vcat(phases["gas"], "unobtainium"))))
        @test startswith(check(model, c3)["error"][1],
                         "ValueError: gas: 'unobtainium' is not registered. Available species: water, hydrogen, ")
        c4 = template_case(temperature_C=[300], phases=merge(phases, Dict("solid" => ["carbon", "Calcite"])))
        @test check(model, c4)["status"] == ["ready"]
        c5 = template_case(temperature_C=[300], carbon_bounds=OrderedDict("gas" => OrderedDict("min" => 0.9),
                                                                         "solid" => OrderedDict("min" => 0.5)))
        @test preflight_messages(model, c5) == ["The carbon shares add up to more than 100 %."]
        c6 = template_case(temperature_C=[300], phases=merge(phases, Dict("gas" => vcat(phases["gas"], "Na+"))))
        t6 = check(model, c6)
        @test t6["error"] == ["ValueError: Ions are only allowed in aqueous: gas: ['Na+']"]
        @test isequal(t6["aqueous_removed"], [missing])   # Python: None
        @test isempty(preflight_messages(model, load_case(SWEEP_TEMPLATE)))
    end

    # One chain at 240 bar through the aqueous cut-off (390 °C): warm starts while
    # the phase selection stays, a fresh LP after it changes.
    chain_case = template_case(temperature_C=[300, 350, 400, 450], pressure_bar=240.0)
    lines = String[]
    warm = run_sweep(model, chain_case; progress=line -> push!(lines, line))

    @testset "continuation: warm starts along a chain" begin
        @test [p[1:2] for p in warm.points] == grid(chain_case)
        @test lines == ["300 °C, 240 bar: solved (start: feasible LP)",
                        "350 °C, 240 bar: solved (start: solution at 300 °C, 240 bar)",
                        "390 °C, 240 bar: solved (start: feasible LP (phase selection changed))",
                        "400 °C, 240 bar: solved (start: solution at 390 °C, 240 bar)",
                        "450 °C, 240 bar: solved (start: solution at 400 °C, 240 bar)"]
        @test all(r -> r isa Result, last.(warm.points))
        @test all(check_solution, results(warm))
        # the warm-started points reach the minimum a cold start finds
        cold = run_sweep(model, replace(chain_case; continuation=false); progress=nothing)
        for (a, b) in zip(results(warm), results(cold))
            @test isapprox(a.gibbs_J_per_mol_C, b.gibbs_J_per_mol_C; rtol=1e-9)
            @test all(isapprox(a.amounts[k], b.amounts[k]; rtol=1e-3, atol=1e-4) for k in keys(b.amounts))
        end
        # Model.sweep is run_sweep
        @test same_results(sweep(model, chain_case; progress=nothing), warm)
        s = summary(warm)
        @test s["status"] == fill("solved", 5) && s["aqueous_removed"] == [false, false, true, true, true]
    end

    @testset "a failed point is recorded and the sweep goes on" begin
        bad = template_case(temperature_C=[300, 400], pressure_bar=240.0, aqueous_cutoff_C=nothing)
        out = run_sweep(model, bad; progress=nothing)
        @test out.points[1][3] isa Result && check_solution(out.points[1][3])
        @test out.points[2][3] isa String && startswith(out.points[2][3], "thermochemistry / aqueous validity: ")
        @test summary(out)["status"] == ["solved", "failed"]
    end

    @testset "threads: a sweep on $(Threads.nthreads()) threads equals the serial one" begin
        # 4 chains of 2 points (warm start at 450 °C), and 8 independent points
        temps = FULL_SWEEP_TESTS ? [300, 350, 450] : [300, 450]
        for c in (template_case(temperature_C=temps, pressure_bar=[240, 260, 280, 300]),
                  template_case(temperature_C=temps, pressure_bar=[240, 260, 280, 300], continuation=false))
            serial = run_sweep(model, c; progress=nothing)
            messages = String[]
            parallel = run_sweep(model, c; progress=line -> push!(messages, line), workers=4)
            @test same_results(parallel, serial)
            @test startswith(messages[1], "Solving $(length(grid(c))) points with 4 workers")
            @test length(messages) == length(grid(c)) + 1
            @test all(check_solution, results(parallel))
        end
        Threads.nthreads() == 1 && @info "sweep threads test ran on 1 thread; set JULIA_NUM_THREADS=4"
    end
end
