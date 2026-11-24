# Optimizer backends, multistart and optimality diagnostics for the Gibbs
# minimization (Python: htg/optimizer.py).
#
# The phase models give G and its exact gradient: for every phase G = sum n_i mu_i
# with mu_i from the same Helmholtz model, so dG/dn_i = mu_i (no finite
# differences). IPOPT also receives the phase-block Hessian, exact (ForwardDiff,
# see phases.jl; Python used forward differences).
#
#     method        "ipopt" (default; Ipopt.jl), "slsqp" (NLopt's LD_SLSQP)
#                   ("trust-constr" and "slsqp-fd" are not available in Julia)
#     starts        total starting points N (default 1): the given start, then the
#                   phase-rich and phase-lean LP starts, then random feasible mixes
#                   (the old `multistart = k` is N = 1 + 2 F + k, still accepted)
#     distinct_composition_tol, distinct_gibbs_tol
#                   converged starts within both limits are the same result; the
#                   distinct results are labelled A, B, ... by increasing G (A is the answer)
#     stability_test tangent-plane test of absent fluid phases (default: on when N > 1)
#
# Thread safety: nothing global is changed. One GibbsProblem (it counts function
# calls and keeps EOS messages) and one Ipopt/NLopt problem per task. The LP
# solver (HiGHS) is called under a lock (HIGHS_LOCK).

"""Python's plain `ValueError`s of optimizer.py and equilibrium.py (same messages)."""
struct EquilibriumError <: ValueError
    msg::String
end

"""Python's `RuntimeError` (a failed solve or check; same messages)."""
struct RuntimeError <: Exception
    msg::String
end
Base.showerror(io::IO, e::RuntimeError) = print(io, e.msg)

"""Solvers of the Python app that the Julia version does not have."""
const UNAVAILABLE_SOLVERS = ("trust-constr", "slsqp-fd")

unavailable_solver_message(method) =
    "The solver $(py_repr(method)) is not available in the Julia version; use 'ipopt' (default) or 'slsqp'."

"""
    GibbsProblem(; phases, pairs, free, phase_indices, ions, T, A_eq, b_eq, ineq_rows,
                 ineq_targets, lower, x_start, G0, scale, model_total_G)

min_x (G(expand(x)) - G0)/scale  s.t.  A x = b,  C x >= d,  x >= lower.

`x` holds the free (phase, species) amounts; `free` and `phase_indices` are
1-based positions in `pairs`. `ineq_rows` (C, one row per inequality) and
`ineq_targets` (d) are `nothing` when there are no inequalities. `phases`
maps phase name => phase model (anything with `gibbs_and_mu`, `hessian`,
`mu_at`, `species`); "solid" is not a fluid phase.
"""
mutable struct GibbsProblem
    phases::OrderedDict{String,Any}
    pairs::Vector{Pair2}
    free::Vector{Int}
    phase_indices::OrderedDict{String,Vector{Int}}
    ions::Set{String}
    T::Float64
    A_eq::Matrix{Float64}
    b_eq::Vector{Float64}
    ineq_rows::Union{Nothing,Matrix{Float64}}
    ineq_targets::Union{Nothing,Vector{Float64}}
    lower::Vector{Float64}
    x_start::Vector{Float64}
    G0::Float64
    scale::Float64
    model_total_G::Any
    eos_errors::Vector{String}
    nfev::Int
end

function GibbsProblem(; phases, pairs, free, phase_indices, ions, T, A_eq, b_eq, ineq_rows=nothing,
                      ineq_targets=nothing, lower, x_start, G0, scale, model_total_G)
    if ineq_rows === nothing || isempty(ineq_rows)
        rows, targets = nothing, nothing
    else
        rows = ineq_rows isa AbstractVector ? reshape(Float64.(ineq_rows), 1, :) : Matrix{Float64}(ineq_rows)
        targets = collect(Float64, ineq_targets isa Number ? [ineq_targets] : ineq_targets)
    end
    A = A_eq isa AbstractVector ? reshape(Float64.(A_eq), 1, :) : Matrix{Float64}(A_eq)
    return GibbsProblem(OrderedDict{String,Any}(String(k) => v for (k, v) in phases),
                        Pair2[(String(p), String(s)) for (p, s) in pairs], collect(Int, free),
                        OrderedDict{String,Vector{Int}}(String(k) => collect(Int, v) for (k, v) in phase_indices),
                        Set{String}(ions), Float64(T), A, collect(Float64, b_eq), rows, targets,
                        collect(Float64, lower), collect(Float64, x_start), Float64(G0), Float64(scale),
                        model_total_G, String[], 0)
end

fluid_phases(problem::GibbsProblem) = [phase for phase in keys(problem.phases) if phase != "solid"]

function expand(problem::GibbsProblem, x)
    amounts = zeros(length(problem.pairs))
    amounts[problem.free] .= x
    return amounts
end

_free_positions(problem::GibbsProblem) = Dict(j => k for (k, j) in enumerate(problem.free))

"""Total G [J/mol feed C] and mu for every (phase, species) pair."""
function gibbs_and_mu(problem::GibbsProblem, amounts)
    amounts = max.(amounts, 0.0)
    G, mu = 0.0, zeros(length(problem.pairs))
    for (phase, model) in problem.phases
        idx = problem.phase_indices[phase]
        g, m = gibbs_and_mu(model, amounts[idx])
        G += g
        mu[idx] .= m
    end
    (!isfinite(G) || any(!isfinite, mu)) && throw(EOSFailure("Nonfinite phase Gibbs energy."))
    return G, mu
end

"""
    fun_and_grad(problem, x; penalty=true) -> (f, grad)

Scaled objective and its exact gradient. An EOS failure returns a large penalty
(penalty=true, for SLSQP whose line search then backtracks) or NaN (for IPOPT,
which cuts the step back).
"""
function fun_and_grad(problem::GibbsProblem, x; penalty=true)
    problem.nfev += 1
    local G, mu
    try
        G, mu = gibbs_and_mu(problem, expand(problem, x))
    catch exc
        exc isa EOSFailure || rethrow()
        length(problem.eos_errors) < 10 && push!(problem.eos_errors, exc.msg)
        penalty || return NaN, fill(NaN, length(x))
        d = x .- problem.x_start
        return 1e6 + sum(abs2, d), 2.0 .* d
    end
    return (G - problem.G0) / problem.scale, mu[problem.free] ./ problem.scale
end

"""Scaled Hessian of the objective as [phase => (free positions, dense block)] (fluid phases only)."""
function hessian_blocks(problem::GibbsProblem, x)
    amounts = max.(expand(problem, x), 0.0)
    position = _free_positions(problem)
    blocks = OrderedDict{String,Tuple{Vector{Int},Matrix{Float64}}}()
    for phase in fluid_phases(problem)
        idx = problem.phase_indices[phase]
        keep = [i for (i, j) in enumerate(idx) if haskey(position, j)]
        isempty(keep) && continue
        H = hessian(problem.phases[phase], amounts[idx])[keep, keep] ./ problem.scale
        blocks[phase] = ([position[idx[i]] for i in keep], H)
    end
    return blocks
end

"""Throw if the optimizer's G (gibbs_and_mu) differs from the model's G (gibbs)."""
function check_against_model(problem::GibbsProblem, x; rtol=1e-9)
    amounts = expand(problem, x)
    ours, _ = gibbs_and_mu(problem, amounts)
    theirs = problem.model_total_G(amounts)
    if abs(ours - theirs) > rtol * max(1.0, abs(theirs))
        throw(RuntimeError("Optimizer G=$(py_repr(ours)) differs from model G=$(py_repr(theirs))"))
    end
    return nothing
end

"""
    kkt(problem, x; major_tol=1e-6) -> OrderedDict

First-order optimality certificate with reduced costs r = (mu - A' lambda)/RT.
At a local minimum r_i = 0 where n_i > lower and r_i >= 0 at the bound.

    descent_RT            max(-r_i): positive means a direction that lowers G
    complementarity_J     RT sum (n_i - lower_i) max(r_i, 0) over minor species
    stationarity_major_RT max |r_i| over species above major_tol
"""
function kkt(problem::GibbsProblem, x; major_tol=1e-6)
    _, mu = gibbs_and_mu(problem, expand(problem, x))
    g = mu[problem.free]
    RT = R * problem.T
    A = problem.A_eq
    # Interior-point solutions end slightly inside an active bound, hence 1e-6.
    if problem.ineq_rows !== nothing
        active = abs.(problem.ineq_rows * x .- problem.ineq_targets) .< 1e-6
        any(active) && (A = vcat(A, problem.ineq_rows[active, :]))
    end
    slack = x .- problem.lower
    major = slack .> major_tol
    fit = count(major) >= size(A, 1) ? major : slack .> 0
    # numpy.linalg.lstsq (rcond=None): minimum-norm least squares
    M = A[:, fit]'
    lam = pinv(M; rtol=eps(Float64) * maximum(size(M); init=1)) * g[fit]
    reduced = (g .- A' * lam) ./ RT
    return OrderedDict{String,Any}(
        "element_potentials_J" => lam[1:size(problem.A_eq, 1)], "reduced_cost_RT" => reduced,
        "tangent_plane_J" => vec(lam' * A),
        "descent_RT" => max(0.0, -minimum(reduced)),
        "complementarity_J" => RT * sum((slack .* max.(reduced, 0))[.!major]; init=0.0),
        "stationarity_major_RT" => any(major) ? maximum(abs.(reduced[major])) : 0.0)
end

# ---------------------------------------------------------------- results of a backend
"""What a backend returns (Python: scipy's OptimizeResult)."""
Base.@kwdef mutable struct OptimizeResult
    x::Vector{Float64}
    fun::Float64 = NaN
    status::Int = 0
    success::Bool = false
    message::String = ""
    nit::Int = 0
end

# ---------------------------------------------------------------- linear programs (HiGHS)
"""HiGHS is called one LP at a time (its scheduler is process-wide); LPs are cheap."""
const HIGHS_LOCK = ReentrantLock()

const _HIGHS_MODEL_STATUS = Dict(
    0 => "Not Set", 1 => "Load error", 2 => "Model error", 3 => "Presolve error", 4 => "Solve error",
    5 => "Postsolve error", 6 => "Empty", 7 => "Optimal", 8 => "Infeasible",
    9 => "Primal infeasible or unbounded", 10 => "Unbounded", 11 => "Bound on objective reached",
    12 => "Target for objective reached", 13 => "Time limit reached", 14 => "Iteration limit reached",
    15 => "Unknown", 16 => "Solution limit reached", 17 => "Interrupted by user")
# SciPy's text for the model status (scipy/optimize/_linprog_highs.py)
const _SCIPY_LP_MESSAGE = Dict(
    7 => "Optimization terminated successfully. ", 13 => "Time limit reached. ",
    14 => "Iteration limit reached. ", 8 => "The problem is infeasible. ", 10 => "The problem is unbounded. ",
    9 => "The problem is unbounded or infeasible. ")
const _SOLUTION_STATUS = Dict(0 => "None", 1 => "Infeasible", 2 => "Feasible")

"""
    linprog(c; A_eq, b_eq, lower, upper=nothing, A_ub=nothing, b_ub=nothing) -> (success, x, message)

min c'x s.t. A_ub x <= b_ub, A_eq x = b_eq, lower <= x <= upper with HiGHS, set
up as SciPy's `linprog(method='highs')` does (inequality rows first, dual
simplex strategy, presolve on). `upper` entries may be Inf. Returns
`x = nothing` when HiGHS does not report an optimal solution; `message`
follows SciPy's text.
"""
function linprog(c; A_eq, b_eq, lower, upper=nothing, A_ub=nothing, b_ub=nothing)
    n = length(c)
    upper = upper === nothing ? fill(Inf, n) : collect(Float64, upper)
    A_ub = A_ub === nothing ? zeros(0, n) : Matrix{Float64}(A_ub)
    b_ub = b_ub === nothing ? Float64[] : collect(Float64, b_ub)
    A = vcat(A_ub, Matrix{Float64}(A_eq))
    row_lower = vcat(fill(-Inf, length(b_ub)), collect(Float64, b_eq))
    row_upper = vcat(b_ub, collect(Float64, b_eq))
    # column-wise (CSC) storage
    starts, index, value = HiGHS.HighsInt[], HiGHS.HighsInt[], Float64[]
    for j in 1:n
        push!(starts, length(value))
        for i in axes(A, 1)
            if A[i, j] != 0
                push!(index, i - 1)
                push!(value, A[i, j])
            end
        end
    end
    m = size(A, 1)
    x = zeros(n); col_dual = zeros(n); row_value = zeros(m); row_dual = zeros(m)
    status, primal = lock(HIGHS_LOCK) do
        highs = HiGHS.Highs_create()
        try
            HiGHS.Highs_setBoolOptionValue(highs, "output_flag", false)
            HiGHS.Highs_setStringOptionValue(highs, "presolve", "on")
            HiGHS.Highs_setIntOptionValue(highs, "simplex_strategy", 1)   # dual, as SciPy
            HiGHS.Highs_passLp(highs, n, m, length(value), HiGHS.kHighsMatrixFormatColwise,
                               HiGHS.kHighsObjSenseMinimize, 0.0, collect(Float64, c), collect(Float64, lower),
                               upper, row_lower, row_upper, starts, index, value)
            HiGHS.Highs_run(highs)
            status = Int(HiGHS.Highs_getModelStatus(highs))
            primal = Ref{HiGHS.HighsInt}(0)
            HiGHS.Highs_getIntInfoValue(highs, "primal_solution_status", primal)
            status == 7 && HiGHS.Highs_getSolution(highs, x, col_dual, row_value, row_dual)
            (status, Int(primal[]))
        finally
            HiGHS.Highs_destroy(highs)
        end
    end
    message = get(_SCIPY_LP_MESSAGE, status, "") * "(HiGHS Status $status: " *
              (status == 7 ? get(_HIGHS_MODEL_STATUS, status, "") :
               "model_status is $(get(_HIGHS_MODEL_STATUS, status, "")); primal_status is " *
               get(_SOLUTION_STATUS, primal, string(primal))) * ")"
    return status == 7, status == 7 ? x : nothing, message
end

# ---------------------------------------------------------------- backends
"""SLSQP of NLopt (`LD_SLSQP`) with the exact gradient (Python: scipy's SLSQP)."""
function _slsqp(problem::GibbsProblem, x0, options)
    n = length(x0)
    opt = NLopt.Opt(:LD_SLSQP, n)
    NLopt.lower_bounds!(opt, problem.lower)
    # SciPy: maxiter iterations; NLopt counts function evaluations.
    NLopt.maxeval!(opt, Int(get(options, "max_iterations", 5000)))
    NLopt.ftol_abs!(opt, Float64(get(options, "ftol", 1e-12)))
    NLopt.min_objective!(opt, function (x, grad)
        f, g = fun_and_grad(problem, x; penalty=true)
        length(grad) > 0 && (grad .= g)
        return f
    end)
    A, b = problem.A_eq, problem.b_eq
    NLopt.equality_constraint!(opt, function (result, x, jac)
        result .= A * x .- b
        length(jac) > 0 && (jac .= A')
        return nothing
    end, fill(1e-12, length(b)))
    if problem.ineq_rows !== nothing
        C, d = problem.ineq_rows, problem.ineq_targets
        # NLopt: c(x) <= 0, here d - C x <= 0
        NLopt.inequality_constraint!(opt, function (result, x, jac)
            result .= d .- C * x
            length(jac) > 0 && (jac .= -C')
            return nothing
        end, fill(1e-12, length(d)))
    end
    start = max.(collect(Float64, x0), problem.lower)   # SciPy clips x0 to the bounds
    minf, minx, ret = NLopt.optimize(opt, start)
    success = ret in (:SUCCESS, :STOPVAL_REACHED, :FTOL_REACHED, :XTOL_REACHED)
    message = success ? "Optimization terminated successfully" :
              ret == :MAXEVAL_REACHED ? "Iteration limit reached" : "NLopt SLSQP stopped: $ret"
    return OptimizeResult(x=collect(Float64, minx), fun=minf, status=success ? 0 : 1, success=success,
                          message=message, nit=NLopt.numevals(opt))
end

const _IPOPT_STATUS = Dict(
    0 => "Solve_Succeeded", 1 => "Solved_To_Acceptable_Level", 2 => "Infeasible_Problem_Detected",
    3 => "Search_Direction_Becomes_Too_Small", 4 => "Diverging_Iterates", 5 => "User_Requested_Stop",
    6 => "Feasible_Point_Found", -1 => "Maximum_Iterations_Exceeded", -2 => "Restoration_Failed",
    -3 => "Error_In_Step_Computation", -4 => "Maximum_CpuTime_Exceeded",
    -10 => "Not_Enough_Degrees_Of_Freedom", -11 => "Invalid_Problem_Definition", -12 => "Invalid_Option",
    -13 => "Invalid_Number_Detected")

"""IPOPT options of the Python app (`_ipopt`)."""
function ipopt_options(options)
    opts = OrderedDict{String,Any}(
        "tol" => 1e-10, "acceptable_tol" => 1e-7, "max_iter" => Int(get(options, "max_iterations", 3000)),
        "print_level" => Int(get(options, "print_level", 0)),
        # The objective is O(1) and starts feasible. A large or adaptively raised
        # barrier parameter pushes trace species into EOS-invalid compositions
        # (observed at 100-200 C with mu_strategy='adaptive'). mu_init = 1e-6
        # sometimes stopped at non-stationary points (manure with ions and
        # minerals, 100-360 C); 1e-9 reached the lowest G found by any setting
        # at every point of the reference cases and example sweeps.
        "mu_strategy" => "monotone", "mu_init" => 1e-9,
        "bound_relax_factor" => 0.0, "honor_original_bounds" => "yes",
        # Default bound_push (1e-2) would move trace species far from an LP vertex.
        "bound_push" => 1e-6, "bound_frac" => 1e-6,
        "jac_c_constant" => "yes", "jac_d_constant" => "yes", "sb" => "yes")
    for (k, v) in get(options, "ipopt_options", Dict())
        opts[String(k)] = v
    end
    return opts
end

"""
IPOPT with the exact gradient and phase-block Hessian (Ipopt.jl).

An EOS failure at a trial point returns NaN, which makes IPOPT cut the step
back. Any other error in a callback stops IPOPT and is thrown again here (as
in Python, where it reaches the caller through juliacall).
"""
function _ipopt(problem::GibbsProblem, x0, options)
    x0 = collect(Float64, x0)
    n = length(x0)
    A = problem.A_eq
    gL, gU = copy(problem.b_eq), copy(problem.b_eq)
    if problem.ineq_rows !== nothing
        A = vcat(A, problem.ineq_rows)
        append!(gL, problem.ineq_targets)
        append!(gU, fill(2e19, length(problem.ineq_targets)))
    end
    Ar, Ac, Av = Int32[], Int32[], Float64[]
    for i in axes(A, 1), j in axes(A, 2)   # numpy.nonzero order (row-major)
        if A[i, j] != 0
            push!(Ar, i); push!(Ac, j); push!(Av, A[i, j])
        end
    end
    # Lower-triangle sparsity of the phase-block Hessian (structure fixed for the solve).
    blocks = hessian_blocks(problem, x0)
    Hr, Hc = Int32[], Int32[]
    for (pos, _) in values(blocks)
        for a in eachindex(pos), b in 1:a
            push!(Hr, max(pos[a], pos[b])); push!(Hc, min(pos[a], pos[b]))
        end
    end
    nit = Ref(0)
    failure = Ref{Any}(nothing)        # an error other than EOSFailure in a callback
    last_x = Ref{Union{Nothing,Vector{Float64}}}(nothing)
    last_fg = Ref{Tuple{Float64,Vector{Float64}}}((NaN, Float64[]))
    function evaluate(x)
        if last_x[] === nothing || last_x[] != x
            last_x[] = copy(x)
            last_fg[] = try
                fun_and_grad(problem, x; penalty=false)
            catch exc
                exc isa InterruptException && rethrow()
                failure[] === nothing && (failure[] = exc)
                (NaN, fill(NaN, length(x)))
            end
        end
        return last_fg[]
    end
    function hess!(out, x)
        nit[] += 1
        k = 0
        try
            for (pos, H) in values(hessian_blocks(problem, x))
                for a in eachindex(pos), b in 1:a
                    out[k += 1] = H[a, b]
                end
            end
        catch exc
            exc isa InterruptException && rethrow()
            exc isa EOSFailure || (failure[] === nothing && (failure[] = exc))
            fill!(out, NaN)
        end
        return out
    end
    eval_f(x) = evaluate(x)[1]
    eval_grad_f(x, g) = (g .= evaluate(x)[2]; nothing)
    function eval_g(x, g)
        fill!(g, 0.0)
        for k in eachindex(Av)
            g[Ar[k]] += Av[k] * x[Ac[k]]
        end
        return nothing
    end
    function eval_jac_g(x, rows, cols, values)
        if values === nothing
            rows .= Ar; cols .= Ac
        else
            values .= Av
        end
        return nothing
    end
    function eval_h(x, rows, cols, sigma, lambda, values)
        if values === nothing
            rows .= Hr; cols .= Hc
        else
            hess!(values, x)
            values .*= sigma
        end
        return nothing
    end
    prob = Ipopt.CreateIpoptProblem(n, problem.lower, fill(2e19, n), length(gL), gL, gU, length(Av), length(Hr),
                                    eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h)
    for (name, value) in ipopt_options(options)
        if value isa AbstractString
            Ipopt.AddIpoptStrOption(prob, name, String(value))
        elseif value isa Integer
            Ipopt.AddIpoptIntOption(prob, name, Int(value))
        else
            Ipopt.AddIpoptNumOption(prob, name, Float64(value))
        end
    end
    # stop at the end of the iteration in which a callback failed
    Ipopt.SetIntermediateCallback(prob, (args...) -> failure[] === nothing)
    prob.x .= x0
    status = Int(Ipopt.IpoptSolve(prob))
    failure[] === nothing || throw(failure[])
    return OptimizeResult(x=copy(prob.x), fun=prob.obj_val, status=status, success=status in (0, 1),
                          message=get(_IPOPT_STATUS, status, string(status)), nit=nit[])
end

"""Optimizer backends by name. Python also had 'trust-constr' (not ported, see UNAVAILABLE_SOLVERS)."""
const BACKENDS = Dict{String,Any}("slsqp" => _slsqp, "ipopt" => _ipopt)

# ---------------------------------------------------------------- multistart
function _lp_point(problem::GibbsProblem, cost)
    A_ub = problem.ineq_rows === nothing ? nothing : -problem.ineq_rows
    b_ub = problem.ineq_rows === nothing ? nothing : -problem.ineq_targets
    success, x, _ = linprog(cost; A_eq=problem.A_eq, b_eq=problem.b_eq, lower=problem.lower, A_ub, b_ub)
    return success ? x : nothing
end

"""The given start, then a phase-rich and a phase-lean LP start for each fluid phase."""
function _phase_starts(problem::GibbsProblem)
    starts = Tuple{String,Vector{Float64}}[("initial", copy(problem.x_start))]
    position = _free_positions(problem)
    for phase in fluid_phases(problem)
        cols = [position[j] for j in problem.phase_indices[phase] if haskey(position, j)]
        isempty(cols) && continue
        for (sign, label) in ((-1.0, "max"), (1.0, "min"))
            cost = zeros(length(problem.free))
            cost[cols] .= sign
            x = _lp_point(problem, cost)
            # Blend toward the given start so no phase is exactly empty.
            x === nothing || push!(starts, ("$label $phase", 0.9 .* x .+ 0.1 .* problem.x_start))
        end
    end
    return starts
end

"""
n_random random mixes of random LP vertices. Julia's RNG (`Xoshiro(seed)`), so
the points differ from Python's NumPy stream; they are deterministic for a seed.
"""
function _random_starts(problem::GibbsProblem, n_random, seed)
    rng = Random.Xoshiro(seed)
    vertices = Vector{Float64}[]
    for _ in 1:max(3, n_random)
        x = _lp_point(problem, randn(rng, length(problem.free)))
        x === nothing || push!(vertices, x)
    end
    isempty(vertices) && return Tuple{String,Vector{Float64}}[]
    starts = Tuple{String,Vector{Float64}}[]
    for k in 0:n_random-1
        e = Random.randexp(rng, length(vertices))      # Dirichlet(1, ..., 1) weights
        w = e ./ sum(e)
        push!(starts, ("random $k", sum(w[i] .* vertices[i] for i in eachindex(vertices))))
    end
    return starts
end

"""
    starting_points(problem; n_random=8, seed=0)

Feasible starts: the given start, phase-rich/phase-lean LP vertices, random
vertex mixes. Every point is a convex combination of LP solutions, so it
satisfies the constraints exactly (to LP precision).
"""
starting_points(problem::GibbsProblem; n_random=8, seed=0) =
    vcat(_phase_starts(problem), _random_starts(problem, n_random, seed))

"""
    candidate_starts(problem, n_total; seed=0)

The first `n_total` starts of the list `starting_points` builds, padded with
random mixes when there are fewer deterministic starts than `n_total`.
n_total = 1 gives only the given start.
"""
function candidate_starts(problem::GibbsProblem, n_total; seed=0)
    n_total = py_int(n_total)
    n_total < 1 && throw(EquilibriumError("starts must be at least 1."))
    n_total == 1 && return Tuple{String,Vector{Float64}}[("initial", copy(problem.x_start))]
    starts = _phase_starts(problem)
    n_total <= length(starts) && return starts[1:n_total]
    return vcat(starts, _random_starts(problem, n_total - length(starts), seed))
end

"""0 -> "A", 25 -> "Z", 26 -> "AA", 27 -> "AB", ..."""
function result_label(k)
    label = ""
    k = Int(k) + 1
    while k > 0
        k, rest = divrem(k - 1, 26)
        label = string(Char('A' + rest)) * label
    end
    return label
end

"""
    group_solutions(solutions; composition_tol=1e-4, gibbs_tol=1.0)

Group converged starts into distinct results.

`solutions` is [(start label, x, G)] with finite G. Two solutions are the same
result when max |dx| <= composition_tol and |dG| <= gibbs_tol. Solutions are
taken by increasing G; each joins the first group whose lowest-G member (its
representative) is within both limits, or starts a new group. Groups are
labelled A, B, ... by increasing G. Returns [OrderedDict("label", "members"
(1-based indices into solutions, by increasing G), "G", "x")].
"""
function group_solutions(solutions; composition_tol=1e-4, gibbs_tol=1.0)
    order = sortperm([Float64(s[3]) for s in solutions])   # stable, as Python's sorted
    groups = OrderedDict{String,Any}[]
    for i in order
        _, x, G = solutions[i]
        x = collect(Float64, x)
        found = false
        for group in groups
            if abs(G - group["G"]) <= gibbs_tol && maximum(abs.(x .- group["x"]); init=0.0) <= composition_tol
                push!(group["members"], i)
                found = true
                break
            end
        end
        found || push!(groups, OrderedDict{String,Any}("members" => [i], "G" => Float64(G), "x" => x))
    end
    return [OrderedDict{String,Any}("label" => result_label(k - 1), "members" => g["members"], "G" => g["G"],
                                    "x" => g["x"]) for (k, g) in enumerate(groups)]
end

# ---------------------------------------------------------------- stability
"""
    tangent_plane_test(problem, x, tangent_plane) -> OrderedDict{phase => min tpd}

Minimum tangent-plane distance (units of RT) for each nearly empty fluid phase.

tpd(y) = sum_i y_i (mu_i(y) - (A' lambda)_i) / RT over compositions y of that
phase. A negative minimum means adding the phase would lower G, so the solution
is not the global minimum of the selected-phase problem. Ions are excluded from
the aqueous trial phase (keeps trial compositions electroneutral). Each search
is NLopt's SLSQP from the uniform composition and from each nearly pure one.
"""
function tangent_plane_test(problem::GibbsProblem, x, tangent_plane)
    amounts = expand(problem, x)
    RT = R * problem.T
    position = _free_positions(problem)
    target = fill(NaN, length(problem.pairs))
    target[problem.free] .= tangent_plane
    report = OrderedDict{String,Float64}()
    for phase in fluid_phases(problem)
        model = problem.phases[phase]
        names = model.species
        idx = problem.phase_indices[phase]
        sum(amounts[idx]) > 1e-6 && continue
        trial = [i for (i, j) in enumerate(idx) if haskey(position, j) && !(names[i] in problem.ions)]
        isempty(trial) && continue
        electrolyte = model isa ElectrolytePhase
        water = findfirst(==("water"), names)
        function tpd(y_trial)
            n = zeros(length(idx))
            n[trial] .= max.(y_trial, 1e-14)
            n ./= sum(n)
            local mu
            try
                if electrolyte && (water === nothing || n[water] <= 1e-9)
                    return Inf, zeros(length(trial))   # the model needs solvent water
                end
                mu = mu_at(model, n)
            catch exc
                exc isa Union{ValueError,ArgumentError} || rethrow()
                return Inf, zeros(length(trial))
            end
            d = (mu .- target[idx])[trial] ./ RT
            return sum(n[trial] .* d), d
        end
        best = Inf
        k = length(trial)
        guesses = vcat([fill(1.0 / k, k)],
                       [[i == j ? 1.0 - 1e-6 * (k - 1) : 1e-6 for i in 1:k] for j in 1:k])
        for y0 in guesses
            isfinite(tpd(y0)[1]) || continue
            opt = NLopt.Opt(:LD_SLSQP, k)
            NLopt.lower_bounds!(opt, zeros(k))
            NLopt.upper_bounds!(opt, ones(k))
            NLopt.maxeval!(opt, 200)
            NLopt.ftol_abs!(opt, 1e-12)
            NLopt.min_objective!(opt, function (y, grad)
                f, d = tpd(y)
                length(grad) > 0 && (grad .= d)
                return f
            end)
            NLopt.equality_constraint!(opt, function (y, grad)
                length(grad) > 0 && (grad .= 1.0)
                return sum(y) - 1.0
            end, 1e-12)
            _, y, _ = NLopt.optimize(opt, y0)
            best = min(best, tpd(y)[1])
        end
        report[phase] = best
    end
    return report
end

const POLISH_THRESHOLD_RT = 0.1

"""
Check IPOPT's answer; polish it with SLSQP when it is not stationary.

IPOPT can report success at a point where species with clearly positive amounts
still have reduced costs of several RT, i.e. not a local minimum (observed for
the ions-and-minerals manure case at 250 C with mu_init = 1e-6). SLSQP with exact
gradients, started from that point, then reaches a lower G. The lower of the two
results is kept.
"""
function _polish(problem::GibbsProblem, result::OptimizeResult, G, options)
    x = max.(result.x, problem.lower)
    threshold = Float64(get(options, "polish_threshold_RT", POLISH_THRESHOLD_RT))
    stationarity = kkt(problem, x)["stationarity_major_RT"]
    info = OrderedDict{String,Any}("stationarity_major_RT" => stationarity, "threshold_RT" => threshold,
                                   "used" => false)
    stationarity <= threshold && return result, info
    local polished, G_new
    try
        polished = _slsqp(problem, x, options)
        G_new = polished.success && all(isfinite, polished.x) ?
                gibbs_and_mu(problem, expand(problem, max.(polished.x, problem.lower)))[1] : NaN
    catch exc
        exc isa InterruptException && rethrow()
        info["message"] = first("$(py_exception_name(exc)): $(_exception_text(exc))", 200)
        return result, info
    end
    info["G_ipopt"] = Float64(G)
    info["G_polished"] = Float64(G_new)
    info["message"] = polished.message
    if isfinite(G_new) && G_new < G
        info["used"] = true
        polished.nit = result.nit + polished.nit
        return polished, info
    end
    return result, info
end

"""Total starting points from the options (`starts`, or the old `multistart`)."""
_total_starts(options) = haskey(options, "starts") ? max(1, py_int(options["starts"])) : nothing

"""
    optimize(problem, options) -> (result, diagnostics)

Run the selected backend from every start; return (OptimizeResult, diagnostics).

Every converged start is kept and grouped into distinct results (see
`group_solutions`). With IPOPT the best solution of each group is checked and
polished; the groups are then sorted by G again and labelled A, B, ... The
returned result is group A. diagnostics["groups"] lists every group with
"label", "G", "x", "n_starts" and "start_labels".

`options` (a dict with string keys): "method", "max_iterations", "starts" (or
the old "multistart"), "seed", "distinct_composition_tol",
"distinct_gibbs_tol", "stability_test", "fallback", "polish",
"polish_threshold_RT", "ftol" (SLSQP), "print_level" and "ipopt_options".
"""
function optimize(problem::GibbsProblem, options)
    method = String(get(options, "method", "ipopt"))
    method in UNAVAILABLE_SOLVERS && throw(EquilibriumError(unavailable_solver_message(method)))
    haskey(BACKENDS, method) ||
        throw(EquilibriumError("Unknown optimizer $(py_repr(method)); choose from $(py_repr(sort(collect(keys(BACKENDS)))))"))
    backend = BACKENDS[method]
    check_against_model(problem, problem.x_start)
    t0 = time()
    seed = get(options, "seed", 0)
    n_total = _total_starts(options)
    if n_total === nothing   # the old option: multistart = k random starts on top of the LP starts
        n_random = py_int(get(options, "multistart", 0))
        starts = n_random > 0 ? starting_points(problem; n_random, seed) :
                 Tuple{String,Vector{Float64}}[("initial", copy(problem.x_start))]
    else
        starts = candidate_starts(problem, n_total; seed)
    end
    runs = OrderedDict{String,Any}[]
    converged = Tuple{String,OptimizeResult,Float64}[]
    result = nothing
    failure = nothing
    for (label, x0) in starts
        start_label = label
        try
            result = backend(problem, x0, options)
        catch exc   # a start may lie outside the EOS domain
            exc isa InterruptException && rethrow()
            result = nothing
            failure = exc
        end
        if method == "ipopt" && py_truth(get(options, "fallback", true)) && (result === nothing || !result.success)
            # IPOPT can stop with a step-computation error when many trial points
            # fail in the equation of state; SLSQP from the same start often succeeds.
            message = result !== nothing ? result.message : "$(py_exception_name(failure)): $(_exception_text(failure))"
            try
                retry = _slsqp(problem, x0, options)
                if retry.success
                    result, label = retry, "$label (ipopt failed: $(first(message, 80)); slsqp)"
                end
            catch exc
                exc isa InterruptException && rethrow()
            end
        end
        if result === nothing
            length(starts) == 1 && throw(failure)
            push!(runs, OrderedDict{String,Any}("start" => label, "G" => NaN,
                                                "message" => first("$(py_exception_name(failure)): $(_exception_text(failure))", 200)))
            continue
        end
        G = NaN
        if result.success && all(isfinite, result.x)
            try
                G, _ = gibbs_and_mu(problem, expand(problem, max.(result.x, problem.lower)))
            catch exc
                exc isa EOSFailure || rethrow()
            end
        end
        push!(runs, OrderedDict{String,Any}("start" => label, "G" => G, "message" => result.message,
                                            "nit" => result.nit))
        isfinite(G) && push!(converged, (start_label, result, Float64(G)))
    end
    diagnostics = OrderedDict{String,Any}("method" => method, "nfev" => problem.nfev, "seconds" => time() - t0,
                                          "starts" => runs)
    if isempty(converged)
        result === nothing && throw(RuntimeError("All optimizer starts failed: $(py_repr(runs))"))
        # return the last (unsuccessful) result so the caller reports its message
        diagnostics["groups"] = OrderedDict{String,Any}[]
        return result, diagnostics
    end
    groups = group_solutions([(label, max.(res.x, problem.lower), G) for (label, res, G) in converged];
                             composition_tol=Float64(get(options, "distinct_composition_tol", 1e-4)),
                             gibbs_tol=Float64(get(options, "distinct_gibbs_tol", 1.0)))
    for group in groups
        best = converged[group["members"][1]]
        group["result"] = best[2]
        group["n_starts"] = length(group["members"])
        group["start_labels"] = [converged[i][1] for i in group["members"]]
        if method == "ipopt" && py_truth(get(options, "polish", true))
            group["result"], group["polish"] = _polish(problem, best[2], best[3], options)
            group["polish"]["used"] && (group["G"] = Float64(group["polish"]["G_polished"]))
        end
        group["x"] = collect(Float64, group["result"].x)
    end
    # Polishing can lower G: sort again (stable, so ties keep their order) and relabel.
    sort!(groups; by=g -> g["G"])
    for (k, group) in enumerate(groups)
        group["label"] = result_label(k - 1)
    end
    result = groups[1]["result"]
    diagnostics["nfev"] = problem.nfev
    diagnostics["seconds"] = time() - t0
    haskey(groups[1], "polish") && (diagnostics["polish"] = groups[1]["polish"])
    diagnostics["groups"] = [OrderedDict{String,Any}(key => group[key] for key in ("label", "G", "x", "n_starts",
                                                                                     "start_labels"))
                             for group in groups]
    if all(isfinite, result.x)
        x = max.(result.x, problem.lower)
        check_against_model(problem, x)
        merge!(diagnostics, kkt(problem, x))
        several = n_total !== nothing ? n_total > 1 : length(starts) > 1
        if py_truth(get(options, "stability_test", several))
            diagnostics["tangent_plane_RT"] = tangent_plane_test(problem, x, diagnostics["tangent_plane_J"])
        end
    end
    return result, diagnostics
end
