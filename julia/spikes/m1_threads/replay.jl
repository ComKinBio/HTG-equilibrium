# Pure-Julia replay of one equilibrium point as the app solves it (M1 spike).
#
# A line-by-line port of the Julia-facing part of the Python package:
#   htg/phases.py     FluidPhase, ElectrolytePhase, SolidPhase, _hessian
#   htg/optimizer.py  GibbsProblem.fun_and_grad / hessian_blocks, _ipopt + _JULIA_IPOPT
# fed with the inputs dump_problems.py captured from the Python app. Clapeyron
# models are built from the same arguments the app passes through juliacall
# (PCSAFT with userlocations; ePCSAFT(components, charges, ideal, PCSAFT,
# hsdDH(components, ConstRSP(eps), refs), refs)), and every Clapeyron call is
# the app's call: fugacity_coefficient(model, P, T, n; phase = root).
#
# Not replayed (pure Python/SciPy in the app, ported later in M2): the LP
# starting points (taken from the dump), SLSQP fallback/polish, KKT and
# tangent-plane diagnostics, standard-state g0 (taken from the dump).

module Replay

using Clapeyron, Ipopt, LinearAlgebra, JSON

const R = 8.31446261815324
const P0 = 100000.0
const MW_WATER_KG = 0.01801528

struct EOSFailure <: Exception
    msg::String
end
struct ValueError <: Exception
    msg::String
end

# ------------------------------------------------------------------ Clapeyron models
function pcsaft_model(spec)
    components = String.(spec["components"])
    ul = Dict{String,Any}()
    for key in ("Mw", "segment", "sigma", "epsilon")
        ul[key] = Float64.(spec[key])
    end
    for key in ("n_H", "n_e")
        haskey(spec, key) && (ul[key] = Int.(spec[key]))
    end
    for key in ("epsilon_assoc", "bondvol")
        if haskey(spec, key)
            d = Dict{Tuple{Tuple{String,String},Tuple{String,String}},Float64}()
            for (a, b, v) in spec[key]
                d[((String(a[1]), String(a[2])), (String(b[1]), String(b[2])))] = Float64(v)
            end
            ul[key] = d
        end
    end
    if haskey(spec, "k")
        ul["k"] = reduce(vcat, permutedims(Float64.(row)) for row in spec["k"])
    end
    verbose = spec["verbose"]
    return verbose === nothing ? PCSAFT(components; userlocations = ul) :
           PCSAFT(components; userlocations = ul, verbose = Bool(verbose))
end

# ------------------------------------------------------------------ phases
abstract type Phase end

struct SolidPhase <: Phase
    mu0::Vector{Float64}
end

mutable struct FluidPhase{M} <: Phase
    mu0::Vector{Float64}
    model::M
    T::Float64
    P::Float64
    root::String
    pure_mu::Union{Nothing,Vector{Float64}}   # per-phase-object cache, as in Python
end

struct ElectrolytePhase{M} <: Phase
    model::M
    T::Float64
    P::Float64
    nall::Int
    split::Int                 # number of neutrals
    positions::Vector{Int}     # 1-based positions of the selected species
    z::Vector{Float64}
    balanced::Vector{Float64}
    mu0_ion::Vector{Float64}
    mu0_neutral::Vector{Float64}
    lnphi_infinity::Vector{Float64}
end

gibbs_and_mu(p::SolidPhase, n) = (dot(n, p.mu0), p.mu0)

function phis(p::FluidPhase, positive)
    φ = fugacity_coefficient(p.model, p.P, p.T, positive; phase = p.root)
    if any(!isfinite, φ) || any(<=(0), φ)
        throw(EOSFailure("Invalid PC-SAFT fugacities: phase=$(p.root)"))
    end
    return φ
end

function mu_at(p::FluidPhase, n)
    positive = max.(n, 1e-14)
    φ = phis(p, positive)
    s = sum(positive)
    return @. p.mu0 + R * p.T * log(positive / s * p.P / P0 * φ)
end

function gibbs_and_mu(p::FluidPhase, n)
    if sum(n) <= 1e-13
        if p.pure_mu === nothing
            k = length(n)
            p.pure_mu = [mu_at(p, [i == j ? 1.0 : 0.0 for j in 1:k])[i] for i in 1:k]
        end
        return 0.0, copy(p.pure_mu)
    end
    mu = mu_at(p, n)
    return dot(n, mu), mu
end

lnphi(p::FluidPhase, n) = log.(phis(p, max.(n, 1e-14)))

function hessian(p::FluidPhase, n)
    sum(n) <= 1e-13 && return zeros(length(n), length(n))
    return _hessian(p, n, sum(max.(n, 1e-14)))
end

function balanced_probe_weights(z, weights)
    w = copy(weights)
    imbalance = dot(z, w)
    candidates = findall(z .* imbalance .< 0)
    if abs(imbalance) > 1e-14
        isempty(candidates) && throw(ValueError("Both positive and negative ionic carriers are needed"))
        j = candidates[1]
        w[j] += -imbalance / z[j]
    end
    return w
end

function phi(p::ElectrolytePhase, n0)
    n = copy(n0)
    (length(n) != p.nall || any(<(0), n)) && throw(ValueError("Invalid aqueous composition vector"))
    s, z = p.split, p.z
    imbalance = dot(n[s+1:end], z)
    smoothing = hypot(imbalance, 2e-14)
    pos, neg = findall(>(0), z), findall(<(0), z)
    (isempty(pos) || isempty(neg)) && throw(ValueError("Aqueous EOS probe requires positive and negative ionic carriers"))
    n[s+pos[1]] += (smoothing - imbalance) / (2 * z[pos[1]])
    n[s+neg[1]] += (smoothing + imbalance) / (-2 * z[neg[1]])
    n[1:s] .+= 1e-14
    n[s+1:end] .+= p.balanced .* 1e-14
    φ = fugacity_coefficient(p.model, p.P, p.T, n; phase = "liquid")
    if length(φ) != length(n) || any(!isfinite, φ) || any(<=(0), φ)
        throw(EOSFailure("Invalid ePC-SAFT fugacity coefficients"))
    end
    return φ
end

function phi_for_optimizer(p::ElectrolytePhase, full)
    try
        return phi(p, full)
    catch e
        e isa ValueError && throw(EOSFailure(e.msg))
        rethrow()
    end
end

function infinite_dilution_reference(p::ElectrolytePhase, water_index)
    scales = (1e-5, 1e-6, 1e-7)
    logs = Vector{Vector{Float64}}()
    for scale in scales
        n = zeros(p.nall)
        n[water_index] = 10.0
        n[p.split+1:end] .= scale .* p.balanced
        push!(logs, log.(phi(p, n)))
    end
    x = sqrt.(collect(scales))
    L = reduce(hcat, logs)'               # 3 x nall
    # np.polyfit(x, y, 1): column-scaled Vandermonde least squares.
    V = hcat(x, ones(3))
    colscale = sqrt.(vec(sum(V .^ 2; dims = 1)))
    c = (V ./ colscale') \ L
    intercept = vec(c[2, :] ./ colscale[2])
    s = p.split
    estimate = logs[3][s+1:end] .- (logs[2][s+1:end] .- logs[3][s+1:end]) .* x[3] ./ (x[2] - x[3])
    maximum(abs.(intercept[s+1:end] .- estimate)) > 0.02 && throw(ValueError("Ion dilution extrapolation did not converge"))
    return intercept
end

function expand(p::ElectrolytePhase, n)
    full = zeros(p.nall)
    full[p.positions] .= n
    return full
end

function gibbs_and_mu(p::ElectrolytePhase, n)
    full = expand(p, n)
    if full[1] <= 1e-9
        sum(full[2:end]) > 1e-10 && throw(ValueError("Aqueous solutes without water"))
        return 0.0, zeros(length(n))
    end
    φ = phi_for_optimizer(p, full)
    RT, s = R * p.T, p.split
    k = p.nall - s
    q = max.(full, 1e-14)
    x = q ./ sum(q)
    mu_neutral = @. p.mu0_neutral + RT * log(x[1:s] * p.P / P0 * φ[1:s])
    molal = max.(full[end-k+1:end], 1e-14) ./ (MW_WATER_KG * full[1])
    mu_ions = @. p.mu0_ion + RT * (log(molal) + log(x[1]) + log(φ[end-k+1:end]) - p.lnphi_infinity[end-k+1:end])
    G = dot(full[1:s], mu_neutral) + dot(full[end-k+1:end], mu_ions)
    return G, vcat(mu_neutral, mu_ions)[p.positions]
end

lnphi(p::ElectrolytePhase, n) = log.(phi_for_optimizer(p, expand(p, n)))[p.positions]

function hessian(p::ElectrolytePhase, n)
    n[1] <= 1e-9 && return zeros(length(n), length(n))
    total = sum(max.(n, 1e-14)) + 1e-14 * (p.nall - length(n))
    return _hessian(p, n, total)
end

function _hessian(p, n, total)
    RT = R * p.T
    alive = n .> 1e-14
    ideal = Diagonal(ifelse.(alive, 1.0 ./ max.(n, 1e-14), 0.0)) .- (alive * alive') ./ total
    base = lnphi(p, n)
    m = length(n)
    jac = zeros(m, m)
    h = 1e-7 * sum(n)
    for j in 1:m
        step = copy(n)
        step[j] += h
        jac[:, j] = (lnphi(p, step) .- base) ./ h
    end
    return RT .* (ideal .+ 0.5 .* (jac .+ jac'))
end

# ------------------------------------------------------------------ building a point
"""Build the phase objects (and their Clapeyron models) of one dumped point."""
function build_phases(case; shared = nothing)
    phases = Pair{String,Phase}[]
    for name in case["phase_order"]
        rec = case["phases"][name]
        t = rec["type"]
        if t == "solid"
            push!(phases, name => SolidPhase(Float64.(rec["mu0"])))
        elseif t == "fluid"
            model = shared === nothing ? pcsaft_model(rec["pcsaft"]) : shared[name]
            push!(phases, name => FluidPhase(Float64.(rec["mu0"]), model, Float64(rec["T"]),
                                             Float64(rec["P"]), String(rec["root"]), nothing))
        else
            model = shared === nothing ? electrolyte_model(rec) : shared[name]
            z = Float64.(rec["z"])
            balanced = balanced_probe_weights(z, Float64.(rec["probe_weights"]))
            all_species = String.(rec["all_species"])
            p = ElectrolytePhase(model, Float64(rec["T"]), Float64(rec["P"]), length(all_species),
                                 length(rec["neutrals"]), Int.(rec["positions"]) .+ 1, z, balanced,
                                 Float64.(rec["mu0_ion"]), Float64.(rec["mu0_neutral"]),
                                 zeros(length(all_species)))
            p.lnphi_infinity .= infinite_dilution_reference(p, findfirst(==("water"), all_species))
            push!(phases, name => p)
        end
    end
    return phases
end

function electrolyte_model(rec)
    components = String.(rec["components"])
    short_range = pcsaft_model(rec["short_range"])
    refs = String[]
    ion_model = hsdDH(components, ConstRSP(Float64(rec["epsilon_r"])), refs)
    return ePCSAFT(components, Int64.(rec["charges"]), short_range.idealmodel, short_range, ion_model, refs)
end

"""Clapeyron models of one point, to share between tasks (shared-model stress test)."""
function build_models(case)
    d = Dict{String,Any}()
    for name in case["phase_order"]
        rec = case["phases"][name]
        rec["type"] == "fluid" && (d[name] = pcsaft_model(rec["pcsaft"]))
        rec["type"] == "electrolyte" && (d[name] = electrolyte_model(rec))
    end
    return d
end

# ------------------------------------------------------------------ the Gibbs problem
mutable struct GibbsProblem
    phases::Vector{Pair{String,Phase}}
    indices::Vector{Vector{Int}}
    npairs::Int
    free::Vector{Int}
    position::Dict{Int,Int}
    G0::Float64
    scale::Float64
    nfev::Int
    neos::Int
end

function GibbsProblem(case, phases)
    free = Int.(case["free"]) .+ 1
    indices = [Int.(case["phase_indices"][name]) .+ 1 for (name, _) in phases]
    GibbsProblem(phases, indices, length(case["pairs"]), free, Dict(j => k for (k, j) in enumerate(free)),
                 Float64(case["G0"]), Float64(case["scale"]), 0, 0)
end

function expand(pr::GibbsProblem, x)
    a = zeros(pr.npairs)
    a[pr.free] .= x
    return a
end

function gibbs_and_mu(pr::GibbsProblem, amounts)
    amounts = max.(amounts, 0.0)
    G, mu = 0.0, zeros(pr.npairs)
    for ((_, ph), idx) in zip(pr.phases, pr.indices)
        g, m = gibbs_and_mu(ph, amounts[idx])
        G += g
        mu[idx] .= m
    end
    (!isfinite(G) || any(!isfinite, mu)) && throw(EOSFailure("Nonfinite phase Gibbs energy."))
    return G, mu
end

function fun_and_grad(pr::GibbsProblem, x)
    pr.nfev += 1
    try
        G, mu = gibbs_and_mu(pr, expand(pr, x))
        return (G - pr.G0) / pr.scale, mu[pr.free] ./ pr.scale
    catch e
        e isa EOSFailure || rethrow()
        pr.neos += 1
        return NaN, fill(NaN, length(x))
    end
end

function hessian_blocks(pr::GibbsProblem, x)
    amounts = max.(expand(pr, x), 0.0)
    blocks = Tuple{Vector{Int},Matrix{Float64}}[]
    for ((name, ph), idx) in zip(pr.phases, pr.indices)
        ph isa SolidPhase && continue
        keep = [i for (i, j) in enumerate(idx) if haskey(pr.position, j)]
        isempty(keep) && continue
        H = hessian(ph, amounts[idx])[keep, keep] ./ pr.scale
        push!(blocks, ([pr.position[idx[i]] for i in keep], H))
    end
    return blocks
end

# ------------------------------------------------------------------ IPOPT (as _ipopt + _JULIA_IPOPT)
struct SolveResult
    x::Vector{Float64}
    status::Int
    obj::Float64
    hessian_evals::Int
    nfev::Int
    eos_failures::Int
    seconds::Float64
end

function ipopt_solve(case, pr::GibbsProblem, x0::Vector{Float64}; extra_options = Pair{String,Any}[])
    t0 = time()
    n = length(x0)
    A = reduce(vcat, permutedims(Float64.(r)) for r in case["A_eq"])
    gL = Float64.(case["b_eq"])
    gU = copy(gL)
    if case["ineq_rows"] !== nothing
        A = vcat(A, reduce(vcat, permutedims(Float64.(r)) for r in case["ineq_rows"]))
        append!(gL, Float64.(case["ineq_targets"]))
        append!(gU, fill(2e19, length(case["ineq_targets"])))
    end
    # np.nonzero order: row-major
    Ar, Ac, Av = Int32[], Int32[], Float64[]
    for i in axes(A, 1), j in axes(A, 2)
        if A[i, j] != 0
            push!(Ar, i); push!(Ac, j); push!(Av, A[i, j])
        end
    end
    Hr, Hc = Int32[], Int32[]
    for (pos, _) in hessian_blocks(pr, x0)
        for a in eachindex(pos), b in 1:a
            push!(Hr, max(pos[a], pos[b])); push!(Hc, min(pos[a], pos[b]))
        end
    end
    nhess = Ref(0)
    last_x = Ref{Vector{Float64}}(Float64[])
    last_fg = Ref{Tuple{Float64,Vector{Float64}}}((NaN, Float64[]))
    function evaluate(x)
        if isempty(last_x[]) || last_x[] != x      # np.array_equal
            last_x[] = copy(x)
            last_fg[] = fun_and_grad(pr, x)
        end
        return last_fg[]
    end
    function hess(x)
        nhess[] += 1
        try
            vals = Float64[]
            for (pos, H) in hessian_blocks(pr, x)
                for a in eachindex(pos), b in 1:a
                    push!(vals, H[a, b])
                end
            end
            return vals
        catch e
            e isa EOSFailure || rethrow()
            return fill(NaN, length(Hr))
        end
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
            values .= sigma .* hess(x)
        end
        return nothing
    end
    prob = Ipopt.CreateIpoptProblem(n, Float64.(case["lower"]), fill(2e19, n), length(gL), gL, gU,
                                    length(Av), length(Hr), eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h)
    for (name, value, kind) in case["ipopt_options"]
        if kind == "str"
            Ipopt.AddIpoptStrOption(prob, name, String(value))
        elseif kind == "int"
            Ipopt.AddIpoptIntOption(prob, name, Int(value))
        else
            Ipopt.AddIpoptNumOption(prob, name, Float64(value))
        end
    end
    for (name, value) in extra_options
        value isa String ? Ipopt.AddIpoptStrOption(prob, name, value) :
        value isa Integer ? Ipopt.AddIpoptIntOption(prob, name, value) :
        Ipopt.AddIpoptNumOption(prob, name, Float64(value))
    end
    prob.x .= x0
    status = Ipopt.IpoptSolve(prob)
    return SolveResult(copy(prob.x), Int(status), prob.obj_val, nhess[], pr.nfev, pr.neos, time() - t0)
end

"""One task = one (point, start): build models + phases, then IPOPT. Returns SolveResult."""
function solve_task(case, start::Int; shared = nothing, extra_options = Pair{String,Any}[])
    t0 = time()
    phases = build_phases(case; shared)
    pr = GibbsProblem(case, phases)
    x0 = Float64.(case["python_ipopt"][start]["x0"])
    r = ipopt_solve(case, pr, x0; extra_options)
    return SolveResult(r.x, r.status, r.obj, r.hessian_evals, r.nfev, r.eos_failures, time() - t0)
end

load_cases(path) = JSON.parsefile(path; dicttype = Dict{String,Any})["cases"]

"""(case index, start index) pairs for every dumped IPOPT start."""
tasks(cases) = [(i, s) for (i, c) in enumerate(cases) for s in eachindex(c["python_ipopt"])]

end # module
