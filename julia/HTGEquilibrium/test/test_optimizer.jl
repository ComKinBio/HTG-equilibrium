# Optimizer: carbon-bound rows, the aqueous cut-off rule, grouping, starting points
# and the multistart bookkeeping with a stand-in problem and scripted backends,
# the LP wrapper and the backends' error handling. Python: tests/test_backend.py
# (the parts without the GUI or result tables).
using Test
using HTGEquilibrium
using HTGEquilibrium: OrderedDict, Phase, EOSFailure, _ipopt, _slsqp, gibbs, gibbs_and_mu, hessian, mu_at
isdefined(Main, :DataTestSetup) || include("data_setup.jl")
using Main.DataTestSetup: REPO

const TEMPLATE = joinpath(REPO, "src", "htg", "templates", "case_template.toml")
template_case(; changes...) = replace(load_case(TEMPLATE); changes...)

const PAIRS = [("gas", "water"), ("gas", "methane"), ("gas", "carbon dioxide"), ("aqueous", "water"),
               ("aqueous", "methane"), ("organic_liquid", "benzene"), ("solid", "carbon"), ("solid", "Magnesite")]
const CARBON = Float64[0, 1, 1, 0, 1, 6, 1, 1]
carbon_mask(phases...) = CARBON .* Float64[p in phases for (p, _) in PAIRS]

# G = c . n (enough for the bookkeeping of starts and results)
struct LinearPhase <: Phase
    species::Vector{String}
    c::Vector{Float64}
end
LinearPhase(c) = LinearPhase(["a", "b"], Float64.(c))
HTGEquilibrium.gibbs_and_mu(p::LinearPhase, n) = (sum(p.c .* n), copy(p.c))
HTGEquilibrium.gibbs(p::LinearPhase, n) = sum(p.c .* n)
HTGEquilibrium.hessian(p::LinearPhase, n) = zeros(2, 2)
HTGEquilibrium.mu_at(p::LinearPhase, n) = copy(p.c)

# A smooth ideal-mixture phase whose G is not linear: G = sum n_i (c_i + log x_i).
# `fail` decides when a composition counts as outside the EOS domain.
struct IdealPhase <: Phase
    species::Vector{String}
    c::Vector{Float64}
    fail::Function
end
function _ideal_mu(p::IdealPhase, n)
    p.fail(n) && throw(EOSFailure("Invalid PC-SAFT fugacities: test"))
    positive = max.(n, 1e-14)
    return p.c .+ log.(positive ./ sum(positive))
end
HTGEquilibrium.gibbs_and_mu(p::IdealPhase, n) = (mu = _ideal_mu(p, n); (sum(n .* mu), mu))
HTGEquilibrium.gibbs(p::IdealPhase, n) = sum(n .* _ideal_mu(p, n))
HTGEquilibrium.mu_at(p::IdealPhase, n) = _ideal_mu(p, n)
function HTGEquilibrium.hessian(p::IdealPhase, n)
    positive = max.(n, 1e-14)
    return [Float64(i == j) / positive[i] - 1 / sum(positive) for i in eachindex(n), j in eachindex(n)]
end

function stand_in_problem(phase=LinearPhase([0.0, 10.0]); x_start=[0.5, 0.5], T=600.0)
    GibbsProblem(phases=OrderedDict("gas" => phase), pairs=[("gas", "a"), ("gas", "b")], free=[1, 2],
                 phase_indices=OrderedDict("gas" => [1, 2]), ions=Set{String}(), T=T, A_eq=[1.0 1.0],
                 b_eq=[1.0], ineq_rows=nothing, ineq_targets=nothing, lower=[0.0, 0.0], x_start=x_start,
                 G0=5.0, scale=1.0, model_total_G=n -> gibbs(phase, n))
end

"""A backend returning the scripted solutions in turn (`nothing`: the start fails)."""
function scripted(answers)
    calls = Vector{Float64}[]
    function backend(problem, x0, options)
        x = answers[length(calls) + 1]
        push!(calls, x0)
        x === nothing && return OptimizeResult(x=x0, success=false, message="failed", nit=1)
        return OptimizeResult(x=Float64.(x), success=true, message="ok", nit=3)
    end
    return backend, calls
end

function with_backend(f, name, backend)
    BACKENDS[name] = backend
    try
        return f()
    finally
        delete!(BACKENDS, name)
    end
end

@testset "optimizer" begin
    @testset "carbon bound rows" begin
        exact, bounds = carbon_bound_rows(OrderedDict("gas" => OrderedDict("min" => 0.6, "max" => 0.9),
                                                      "solid" => OrderedDict("max" => 0.05),
                                                      "aqueous" => OrderedDict("exact" => 0.03),
                                                      "organic_liquid" => OrderedDict("min" => 0.01)),
                                          PAIRS, CARBON, "char")
        @test [(p, f) for (p, _, f) in exact] == [("aqueous", 0.03)]
        @test exact[1][2] == carbon_mask("aqueous")
        @test [(p, k, f) for (p, k, _, f) in bounds] ==
              [("gas", "min", 0.6), ("gas", "max", 0.9), ("solid", "max", 0.05), ("organic_liquid", "min", 0.01)]
        @test bounds[1][3] == carbon_mask("gas")
        @test bounds[4][3] == carbon_mask("organic_liquid")
        char = zeros(length(PAIRS)); char[findfirst(==(("solid", "carbon")), PAIRS)] = 1.0
        @test bounds[3][3] == char   # char only

        # the solid row follows the target
        _, bounds = carbon_bound_rows(OrderedDict("solid" => OrderedDict("min" => 0.1)), PAIRS, CARBON, "all_solids")
        @test bounds[1][3] == carbon_mask("solid")
        err = try
            carbon_bound_rows(OrderedDict("solid" => OrderedDict("exact" => 0.1)), PAIRS[1:end-2], CARBON[1:end-2], "char")
        catch e
            e
        end
        @test err isa EquilibriumError && err isa ValueError && occursin("Select solid 'carbon'", err.msg)
        @test_throws EquilibriumError carbon_bound_rows(OrderedDict("gas" => OrderedDict("min" => 1.5)), PAIRS,
                                                        CARBON, "char")

        # the old gas minimum gives the old row: gas carbon >= fraction * feed carbon
        old = template_case(gas_carbon_min_fraction=0.9)
        @test old.carbon_bounds == OrderedDict("gas" => OrderedDict("min" => 0.9))
        exact, bounds = carbon_bound_rows(old.carbon_bounds, PAIRS, CARBON, old.solid_carbon_target)
        @test isempty(exact) && length(bounds) == 1
        @test bounds[1][[1, 2, 4]] == ("gas", "min", 0.9)
        @test bounds[1][3] == carbon_mask("gas")
        exact, _ = carbon_bound_rows(template_case(gas_carbon_fraction=0.8, solid_carbon_fraction=0.1,
                                                   phase_carbon_fractions=OrderedDict("organic_liquid" => 0.05)).carbon_bounds,
                                     PAIRS, CARBON, "char")
        @test [(p, f) for (p, _, f) in exact] == [("gas", 0.8), ("solid", 0.1), ("organic_liquid", 0.05)]
        @test carbon_bound_violation([("gas", "min", carbon_mask("gas"), 0.5), ("solid", "exact", char, 0.1)],
                                     [0, 0.2, 0.2, 0, 0, 0, 0.3, 0], 1.0) ≈ 0.2
    end

    @testset "aqueous cut-off drops the aqueous bound" begin
        registry = (ion=Dict{String,Any}(),)   # stand-in: apply_phase_rules reads only the ion table
        c = template_case(carbon_bounds=OrderedDict("aqueous" => OrderedDict("exact" => 0.03),
                                                    "gas" => OrderedDict("min" => 0.5)), aqueous_cutoff_C=390)
        for (T, dropped) in ((350.0, String[]), (390.0, ["aqueous carbon"]), (450.0, ["aqueous carbon"]))
            phases, minimum, _, liquid, gas_water, report = apply_phase_rules(at(c, T), registry)
            @test report["dropped_constraints"] == dropped
            @test report["aqueous_disabled"] == (T >= 390)
            @test isempty(phases["aqueous"]) == (T >= 390)
            @test report["removed_aqueous_species"] == (T >= 390 ? c.phases["aqueous"] : String[])
            T >= 390 && @test !liquid && gas_water && all(k -> k[1] != "aqueous", keys(minimum))
        end
        bad = template_case(phases=OrderedDict("gas" => ["water", "Na+"], "aqueous" => ["water"]))
        err = try apply_phase_rules(at(bad, 300.0), (ion=Dict("Na+" => nothing),)) catch e e end
        @test err isa EquilibriumError && err.msg == "Ions are only allowed in aqueous: gas: ['Na+']"
    end

    @testset "grouping" begin
        x = [1.0, 0.0, 0.5]
        runs = [("initial", x, -100.0), ("max gas", x .+ 5e-5, -100.4), ("random 0", x .+ [0, 0, 1e-3], -100.2),
                ("random 1", x, -98.0), ("min gas", x .+ 2e-5, -100.39)]
        groups = group_solutions(runs; composition_tol=1e-4, gibbs_tol=1.0)
        # lowest G first: 'max gas' (-100.4); 'min gas', 'initial' are within both limits of it;
        # 'random 0' differs in composition; 'random 1' differs in G only
        @test [g["label"] for g in groups] == ["A", "B", "C"]
        @test [[runs[i][1] for i in g["members"]] for g in groups] ==
              [["max gas", "min gas", "initial"], ["random 0"], ["random 1"]]
        @test [g["G"] for g in groups] == [-100.4, -100.2, -98.0]
        loose = group_solutions(runs; composition_tol=1e-2, gibbs_tol=5.0)
        @test length(loose) == 1 && length(loose[1]["members"]) == 5
        tight_g = group_solutions(runs; composition_tol=1e-2, gibbs_tol=0.05)
        @test [[runs[i][1] for i in g["members"]] for g in tight_g] ==
              [["max gas", "min gas"], ["random 0"], ["initial"], ["random 1"]]

        # identical runs collapse
        y = [0.2, 0.3]
        groups = group_solutions([("a", y, -5.0), ("b", copy(y), -5.0), ("c", copy(y), -5.0)])
        @test length(groups) == 1 && groups[1]["members"] == [1, 2, 3] && groups[1]["label"] == "A"
        @test result_label(0) == "A" && result_label(25) == "Z" && result_label(26) == "AA" &&
              result_label(27) == "AB" && result_label(26 + 26 * 26) == "AAA"
    end

    @testset "every start is kept and grouped" begin
        backend, calls = scripted([[0.99995, 0.00005], [1.0, 0.0], nothing, [0.0, 1.0], [0.3, 0.7]])
        result, diagnostics = with_backend("scripted", backend) do
            optimize(stand_in_problem(), Dict("method" => "scripted", "starts" => 5, "seed" => 1))
        end
        @test length(calls) == 5
        labels = [run["start"] for run in diagnostics["starts"]]
        @test labels == ["initial", "max gas", "min gas", "random 0", "random 1"]
        @test isnan(diagnostics["starts"][3]["G"]) && diagnostics["starts"][3]["message"] == "failed"
        groups = diagnostics["groups"]
        @test [g["label"] for g in groups] == ["A", "B", "C"]
        @test groups[1]["start_labels"] == ["max gas", "initial"] && groups[1]["n_starts"] == 2
        @test groups[2]["G"] ≈ 7.0 && groups[3]["G"] ≈ 10.0
        @test result.x ≈ [1.0, 0.0]   # the main result is A
        @test haskey(diagnostics, "tangent_plane_RT")   # stability test on by default with several starts
        @test diagnostics["method"] == "scripted" && haskey(diagnostics, "descent_RT")
    end

    @testset "one start is the old behaviour" begin
        backend, calls = scripted([[1.0, 0.0]])
        result, diagnostics = with_backend(() -> optimize(stand_in_problem(), Dict("method" => "scripted")),
                                           "scripted", backend)
        @test length(calls) == 1 && [run["start"] for run in diagnostics["starts"]] == ["initial"]
        @test [g["label"] for g in diagnostics["groups"]] == ["A"]
        @test !haskey(diagnostics, "tangent_plane_RT")
        # the old option: multistart = k random starts on top of the LP starts
        backend, calls = scripted(fill([1.0, 0.0], 5))
        _, diagnostics = with_backend(() -> optimize(stand_in_problem(), Dict("method" => "scripted", "multistart" => 2)),
                                      "scripted", backend)
        @test [run["start"] for run in diagnostics["starts"]] == ["initial", "max gas", "min gas", "random 0", "random 1"]
        @test haskey(diagnostics, "tangent_plane_RT")
    end

    @testset "candidate starts are padded with random mixes" begin
        problem = stand_in_problem()
        @test [label for (label, _) in candidate_starts(problem, 1)] == ["initial"]
        @test [label for (label, _) in candidate_starts(problem, 2)] == ["initial", "max gas"]
        starts = candidate_starts(problem, 6; seed=3)
        @test [label for (label, _) in starts] == ["initial", "max gas", "min gas", "random 0", "random 1", "random 2"]
        old = starting_points(problem; n_random=3, seed=3)   # the old multistart = 3
        @test [label for (label, _) in old] == [label for (label, _) in starts]
        @test all(a == b for ((_, a), (_, b)) in zip(old, starts))
        for (_, x) in starts
            @test sum(x) ≈ 1.0 && all(x .>= -1e-12)
        end
        # LP starts blend 0.9 of the vertex with 0.1 of the given start
        @test starts[2][2] ≈ [0.95, 0.05] || starts[2][2] ≈ [0.05, 0.95]
        # deterministic for a seed (Julia's RNG; the points differ from NumPy's)
        @test candidate_starts(problem, 6; seed=3)[4:6] == starts[4:6]
        @test candidate_starts(problem, 6; seed=4)[4][2] != starts[4][2]
        @test_throws EquilibriumError candidate_starts(problem, 0)
    end

    @testset "linprog (HiGHS, SciPy's setup)" begin
        ok, x, message = linprog([1.0, 2.0]; A_eq=[1.0 1.0], b_eq=[1.0], lower=[0.0, 0.0])
        @test ok && x ≈ [1.0, 0.0] && startswith(message, "Optimization terminated successfully. (HiGHS Status 7")
        ok, x, _ = linprog([1.0, 0.0]; A_eq=[1.0 1.0], b_eq=[1.0], lower=[0.0, 0.0], A_ub=[-1.0 0.0], b_ub=[-0.25])
        @test ok && x ≈ [0.25, 0.75]
        ok, x, message = linprog([0.0, 0.0]; A_eq=[1.0 1.0], b_eq=[-1.0], lower=[0.0, 0.0])
        @test !ok && x === nothing
        @test startswith(message, "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible")
        # upper bounds
        ok, x, _ = linprog([-1.0, 0.0]; A_eq=[1.0 1.0], b_eq=[1.0], lower=[0.0, 0.0], upper=[0.3, Inf])
        @test ok && x ≈ [0.3, 0.7]
    end

    @testset "IPOPT and SLSQP on a smooth stand-in" begin
        # min sum n_i (c_i + log x_i) with n_a + n_b = 1: x_i = exp(-c_i) / sum exp(-c)
        phase = IdealPhase(["a", "b"], [0.0, 1.0], n -> false)
        expected = exp.(-phase.c) ./ sum(exp.(-phase.c))
        for method in ("ipopt", "slsqp")
            result, d = optimize(stand_in_problem(phase), Dict("method" => method))
            @test result.success
            @test result.x ≈ expected rtol = 1e-7
            @test d["stationarity_major_RT"] < 1e-6
            @test d["element_potentials_J"] ≈ [-log(sum(exp.(-phase.c)))] rtol = 1e-6
            method == "ipopt" && @test !d["polish"]["used"]
        end
        r = _ipopt(stand_in_problem(phase), [0.5, 0.5], Dict())
        @test r.message == "Solve_Succeeded" && r.status == 0 && r.nit > 0
        r = _slsqp(stand_in_problem(phase), [0.5, 0.5], Dict())
        @test r.message == "Optimization terminated successfully"
        @test r.x ≈ expected rtol = 1e-7
        # IPOPT options: the app's defaults, overridden by ipopt_options
        opts = HTGEquilibrium.ipopt_options(Dict("ipopt_options" => Dict("mu_init" => 1e-6), "max_iterations" => 7))
        @test opts["mu_init"] == 1e-6 && opts["max_iter"] == 7 && opts["mu_strategy"] == "monotone" &&
              opts["tol"] == 1e-10
    end

    @testset "EOS failures and errors in the backends" begin
        # EOS failure: NaN for IPOPT, a penalty for SLSQP; messages are kept (at most 10)
        phase = IdealPhase(["a", "b"], [0.0, 1.0], n -> n[1] > 0.9)
        problem = stand_in_problem(phase)
        f, g = fun_and_grad(problem, [0.95, 0.05]; penalty=false)
        @test isnan(f) && all(isnan, g)
        f, g = fun_and_grad(problem, [0.95, 0.05])
        @test f ≈ 1e6 + 2 * 0.45^2 && g ≈ 2 .* ([0.95, 0.05] .- 0.5)
        @test problem.eos_errors == fill("Invalid PC-SAFT fugacities: test", 2) && problem.nfev == 2
        result, _ = optimize(stand_in_problem(phase), Dict("method" => "ipopt"))
        @test result.x ≈ exp.(-phase.c) ./ sum(exp.(-phase.c)) rtol = 1e-7   # the minimum is inside the domain

        # another error inside an IPOPT callback stops IPOPT and is thrown again
        boom = IdealPhase(["a", "b"], [0.0, 1.0], n -> n[1] > 0.6 && error("boom"))
        @test_throws ErrorException("boom") _ipopt(stand_in_problem(boom; x_start=[0.5, 0.5]), [0.5, 0.5], Dict())
        # ... and with one start (and the SLSQP fallback failing too) optimize throws it
        @test_throws ErrorException("boom") optimize(stand_in_problem(boom), Dict("method" => "ipopt"))
        # with several starts the failed start is recorded
        backend(problem, x0, options) = x0 != problem.x_start ? error("bad start") : _ipopt(problem, x0, options)
        _, d = with_backend(() -> optimize(stand_in_problem(phase), Dict("method" => "fragile", "starts" => 3)),
                            "fragile", backend)
        @test d["starts"][2]["start"] == "max gas" && d["starts"][2]["message"] == "ErrorException: bad start"
        @test isnan(d["starts"][2]["G"]) && isfinite(d["starts"][1]["G"])
    end

    @testset "solvers that are not ported" begin
        problem = stand_in_problem()
        for method in ("trust-constr", "slsqp-fd")
            err = try optimize(problem, Dict("method" => method)) catch e e end
            @test err isa EquilibriumError
            @test err.msg == "The solver '$method' is not available in the Julia version; use 'ipopt' (default) or 'slsqp'."
        end
        err = try optimize(problem, Dict("method" => "newton")) catch e e end
        @test err.msg == "Unknown optimizer 'newton'; choose from ['ipopt', 'slsqp']"
        @test UNAVAILABLE_SOLVERS == ("trust-constr", "slsqp-fd")
    end

    @testset "check against the model and the KKT certificate" begin
        phase = IdealPhase(["a", "b"], [0.0, 1.0], n -> false)
        problem = stand_in_problem(phase)
        @test check_against_model(problem, [0.3, 0.7]) === nothing
        problem.model_total_G = n -> gibbs(phase, n) + 1.0
        err = try check_against_model(problem, [0.3, 0.7]) catch e e end
        @test err isa RuntimeError && startswith(err.msg, "Optimizer G=")
        k = kkt(stand_in_problem(phase), [0.5, 0.5])
        @test k["descent_RT"] > 0 && k["stationarity_major_RT"] > 0
        @test k["tangent_plane_J"] ≈ fill(k["element_potentials_J"][1], 2)
        blocks = hessian_blocks(stand_in_problem(phase), [0.25, 0.75])
        @test collect(keys(blocks)) == ["gas"] && blocks["gas"][1] == [1, 2]
        @test blocks["gas"][2] ≈ [1/0.25-1 -1; -1 1/0.75-1]
    end
end
