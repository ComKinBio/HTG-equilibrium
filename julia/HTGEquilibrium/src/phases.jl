# Phase models: Gibbs energy G(n) and chemical potentials mu(n) of one phase
# (Python: htg/phases.py).
#
# Every phase is built for fixed T [K] and P [Pa] and has
#
#     gibbs(phase, n)          G [J] of amounts n [mol], ordered like `phase.species`
#     gibbs_and_mu(phase, n)   (G, mu): mu_i = dG/dn_i; G = sum n_i mu_i (Euler)
#     hessian(phase, n)        dmu/dn (fluid and electrolyte phases)
#
# Formulas
#     Fluid molecules (PC-SAFT):  mu_i = g°_i(T) + RT ln(x_i P phi_i / P0)
#     Pure solids:                mu_i = g°_i(T, P)
#     Aqueous electrolyte (ePC-SAFT + HKF):
#         neutrals  mu_i = g°_i(T) + RT ln(x_i P phi_i / P0)
#         ions      mu_i = g°_i,HKF(T, P) + RT [ln m_i + ln x_w + ln phi_i - ln phi_i^inf]
#     with molality m_i = n_i / (M_w n_w) and phi_i^inf the infinite-dilution limit.
#
# Clapeyron is used directly: the models are built with the arguments the
# Python app passes through juliacall (PCSAFT(components; userlocations);
# ePCSAFT(components, charges, ideal, PCSAFT, hsdDH(components, ConstRSP(eps_r), refs), refs))
# and the only call while solving is fugacity_coefficient(model, P, T, n; phase).
#
# Threads: nothing here is global or cached between phase objects. A Clapeyron
# model is only read, so one model can be shared by several phase objects and
# threads (M1). An ElectrolytePhase is not changed after construction; the only
# cache, FluidPhase's pure-component mu (computed once, as in Python), is filled
# under the phase's own lock. The intended pattern is still one set of phase
# objects per task.
#
# The functions take any real element type, so ForwardDiff can differentiate
# through them (used by the tests to check mu and the Hessian).

"""
    EOSFailure(msg)

The equation of state gave no valid fugacity coefficients at this composition
(Python: `class EOSFailure(ValueError)`). The optimizer treats it as a failed
step and cuts the step back.
"""
struct EOSFailure <: ValueError
    msg::String
end

"""A ValueError raised by the phase models (invalid input, state or parameters); Python: ValueError."""
struct PhaseError <: ValueError
    msg::String
end

"""Abstract type of the phase models (`SolidPhase`, `FluidPhase`, `ElectrolytePhase`)."""
abstract type Phase end

_dot(a, b) = sum(a[i] * b[i] for i in eachindex(a, b); init=zero(eltype(a)) * zero(eltype(b)))
_plain(x) = ForwardDiff.value.(x)   # Float64 values for messages (Duals when differentiated)

"""
    association_sites(association_scheme, epsilonAB, kappaAB) -> (n_H, n_e)

Association sites of a molecule: '2B' (1, 1), '3B' (1, 2); otherwise (1, 1) when
both association parameters are finite and positive, else (0, 0).
"""
function association_sites(association_scheme, epsilonAB, kappaAB)
    scheme = uppercase(strip(py_str(association_scheme)))
    scheme == "2B" && return 1, 1
    scheme == "3B" && return 1, 2
    if isfinite(epsilonAB) && isfinite(kappaAB) && epsilonAB > 0.0 && kappaAB > 0.0
        return 1, 1
    end
    return 0, 0
end

function _first_value(mapping, names, default=NaN)
    for name in names
        if haskey(mapping, name)
            value = mapping[name]
            isna(value) || return value
        end
    end
    return default
end

function _required_float(value, parameter, species)
    result = py_tryfloat(value)
    result === nothing && throw(PhaseError("$species: invalid $parameter = $(py_repr(value))"))
    isfinite(result) || throw(PhaseError("$species: missing $parameter"))
    return result
end

_as_float(value) = something(py_tryfloat(value), NaN)

const AssocKey = Tuple{Tuple{String,String},Tuple{String,String}}

"""
    build_pcsaft_model(species, resolved, binary_interactions) -> Clapeyron.PCSAFT

Clapeyron PCSAFT model of the molecules `species` from their resolved
parameter sets (`resolve_parameters`), with the registry's binary interactions
(`registry.binary_interactions`) between the components present.
"""
function build_pcsaft_model(species, resolved, binary_interactions)
    components, MW, segment, sigma, epsilon = String[], Float64[], Float64[], Float64[], Float64[]
    n_H, n_e = Int[], Int[]
    epsilon_assoc, bondvol = Dict{AssocKey,Float64}(), Dict{AssocKey,Float64}()
    for name in species
        haskey(resolved, name) || throw(KeyError("No resolved PC-SAFT parameter set for '$name'."))
        p = resolved[name]
        component = String(strip(py_str(p["Formula"])))
        push!(components, component)
        push!(MW, molecular_weight(component, ATOMIC_WEIGHTS))
        push!(segment, py_float(p["m"]))
        push!(sigma, py_float(p["sigma"]))
        push!(epsilon, py_float(p["epsilon"]))
        epsilonAB = _as_float(get(p, "epsilonAB", NaN))
        kappaAB = _as_float(get(p, "kappaAB", NaN))
        nH, ne = association_sites(get(p, "association_scheme", "-"), epsilonAB, kappaAB)
        push!(n_H, nH)
        push!(n_e, ne)
        if nH > 0 && ne > 0 && isfinite(epsilonAB) && isfinite(kappaAB)
            key = ((component, "e"), (component, "H"))
            epsilon_assoc[key] = epsilonAB
            bondvol[key] = kappaAB
        end
    end
    for (label, values) in ("Mw" => MW, "segment" => segment, "sigma" => sigma, "epsilon" => epsilon)
        all(isfinite, values) || throw(PhaseError("$label contains non-finite values."))
    end
    userlocations = Dict{String,Any}("Mw" => MW, "segment" => segment, "sigma" => sigma, "epsilon" => epsilon)
    if any(>(0), n_H) || any(>(0), n_e)
        userlocations["n_H"] = n_H
        userlocations["n_e"] = n_e
        isempty(epsilon_assoc) || (userlocations["epsilon_assoc"] = epsilon_assoc)
        isempty(bondvol) || (userlocations["bondvol"] = bondvol)
    end
    k_matrix = nothing
    for ((a, b), kij) in binary_interactions
        if a in components && b in components
            k_matrix === nothing && (k_matrix = zeros(length(components), length(components)))
            i, j = findfirst(==(a), components), findfirst(==(b), components)
            k_matrix[i, j] = k_matrix[j, i] = kij
        end
    end
    k_matrix === nothing || (userlocations["k"] = k_matrix)
    return PCSAFT(components; userlocations=userlocations)
end

# ---------------------------------------------------------------- solids
"""
    SolidPhase(species, mu0)

Pure solids: G = sum n_i g°_i(T, P).
"""
struct SolidPhase <: Phase
    species::Vector{String}
    mu0::Vector{Float64}
    SolidPhase(species, mu0) = new(collect(String, species), collect(Float64, mu0))
end

gibbs(p::SolidPhase, n) = _dot(n, p.mu0)
gibbs_and_mu(p::SolidPhase, n) = (_dot(n, p.mu0), copy(p.mu0))

# ---------------------------------------------------------------- fluids
"""
    FluidPhase(name, species, mu0, eos_model, T, P, root)

Molecular fluid phase with PC-SAFT fugacity coefficients. `root`: "vapour"
(gas) or "liquid" (organic liquid, molecular aqueous). The root choice alone
does not establish phase stability. `eos_model` (e.g. from `build_pcsaft_model`)
is only read and may be shared with other phase objects and threads.
"""
struct FluidPhase{M} <: Phase
    name::String
    species::Vector{String}
    mu0::Vector{Float64}
    eos_model::M
    T::Float64
    P::Float64
    root::String
    _pure_mu::Base.RefValue{Union{Nothing,Vector{Float64}}}
    _lock::ReentrantLock
end
FluidPhase(name, species, mu0, eos_model, T, P, root) =
    FluidPhase(String(name), collect(String, species), collect(Float64, mu0), eos_model, Float64(T),
               Float64(P), String(root), Ref{Union{Nothing,Vector{Float64}}}(nothing), ReentrantLock())

function _phis(p::FluidPhase, positive)
    phis = fugacity_coefficient(p.eos_model, p.P, p.T, positive; phase=p.root)
    if any(!isfinite, phis) || any(<=(0), phis)
        throw(EOSFailure("Invalid PC-SAFT fugacities: phase=$(p.root), " *
                         "T=$(py_format_g(p.T - 273.15, 6)) C, P=$(py_format_g(p.P / 100000.0, 6)) bar, " *
                         "species=$(py_repr(p.species)), x=$(py_repr(_plain(positive ./ sum(positive))))"))
    end
    return phis
end

"""Chemical potentials [J/mol] at amounts `n` (no special case for an empty phase)."""
function mu_at(p::FluidPhase, n)
    positive = max.(n, 1e-14)
    phis = _phis(p, positive)
    total = sum(positive)
    return @. p.mu0 + R * p.T * log(positive / total * p.P / P0 * phis)
end

function gibbs(p::FluidPhase, n)
    sum(n) <= 1e-13 && return zero(float(eltype(n)))
    return _dot(n, mu_at(p, n))
end

function gibbs_and_mu(p::FluidPhase, n)
    if sum(n) <= 1e-13
        # G = 0 for an empty phase. G is homogeneous of degree one, so the
        # one-sided derivative along e_i is the pure-i chemical potential.
        pure = lock(p._lock) do
            if p._pure_mu[] === nothing
                k = length(p.species)
                p._pure_mu[] = [mu_at(p, [Float64(i == j) for j in 1:k])[i] for i in 1:k]
            end
            p._pure_mu[]
        end
        return 0.0, copy(pure)
    end
    mu = mu_at(p, n)
    return _dot(n, mu), mu
end

"""ln of the fugacity coefficients at amounts `n`."""
lnphi(p::FluidPhase, n) = log.(_phis(p, max.(n, 1e-14)))

function hessian(p::FluidPhase, n)
    sum(n) <= 1e-13 && return zeros(length(n), length(n))
    return _hessian(p, n, sum(max.(n, 1e-14)))
end

# ---------------------------------------------------------------- aqueous electrolyte
"""
    ElectrolytePhase(water_model, model, selected, T, P, resolved, atomic_weights)

Aqueous electrolyte: ePC-SAFT (PC-SAFT + Debye-Hückel) with HKF ion standard
states. Python: `ElectrolytePhase(model, selected, T, P, resolved, atomic_weights)`;
the Julia version takes the case's water model first (`water_model_of(case)`),
used for the validity check, the dielectric constant and the HKF g° of the
ions, so a phase never mixes two water models.

The equation of state always contains the registry's reference neutrals and
ions (in that order) followed by the selected species; unselected components
carry trace probe amounts. `species` lists the case's selected aqueous
species, `all_species` the full component vector, `positions` the (1-based)
positions of `species` in `all_species`. `resolved` holds the PC-SAFT
parameters of the neutrals (`resolve_parameters`).
"""
struct ElectrolytePhase{W<:WaterModel,S,M} <: Phase
    water_model::W
    T::Float64
    P::Float64
    species::Vector{String}
    neutrals::Vector{String}
    ions::Vector{String}
    all_species::Vector{String}
    positions::Vector{Int}
    z::Vector{Float64}
    probe_weights::OrderedDict{String,Float64}
    water::WaterProperties
    components::Vector{String}
    charges::Vector{Int}
    short_range::S
    eos_model::M
    mu0_ion::Vector{Float64}
    mu0_neutral::Vector{Float64}
    lnphi_infinity::Vector{Float64}
end

function ElectrolytePhase(water_model::WaterModel, model::Model, selected, T, P, resolved, atomic_weights)
    registry = model.registry
    ions_table = registry.ion
    T, P = Float64(T), Float64(P)
    species = collect(String, selected)
    neutrals = unique(vcat(registry.aqueous_reference_neutrals, [s for s in species if !haskey(ions_table, s)]))
    ions = unique(vcat(registry.aqueous_reference_ions, [s for s in species if haskey(ions_table, s)]))
    all_species = vcat(neutrals, ions)
    first(all_species) == "water" || throw(PhaseError("The aqueous reference species must start with water."))
    positions = [findfirst(==(s), all_species) for s in species]
    z = Float64[ions_table[s]["charge"] for s in ions]
    probe_weights = registry.probe_weights

    water = water_properties(water_model, T, P)
    eps = water.epsilon
    if !(isfinite(water.rho) && water.rho >= 350 && isfinite(eps) && eps > 1)
        throw(PhaseError("State outside SUPCRT aqueous validity checks: " *
                         "rho=$(py_repr(water.rho)), epsilon=$(py_repr(eps))"))
    end
    components, charges, short_range = _build_short_range(model, neutrals, ions, all_species, resolved,
                                                          atomic_weights)
    refs = String[]
    ion_model = hsdDH(components, ConstRSP(Float64(eps)), refs)
    eos_model = ePCSAFT(components, charges, short_range.idealmodel, short_range, ion_model, refs)
    mu0_ion = hkf_ion_g0(water_model, model, ions, T, P)
    mu0_neutral = nasa_g0(model, neutrals, T)
    phase = ElectrolytePhase(water_model, T, P, species, neutrals, ions, all_species, positions, z,
                             probe_weights, water, components, charges, short_range, eos_model,
                             mu0_ion, mu0_neutral, zeros(length(all_species)))
    phase.lnphi_infinity .= _infinite_dilution_reference(phase)
    return phase
end

function _build_short_range(model::Model, neutrals, ions, all_species, resolved, atomic_weights)
    registry = model.registry
    molecular = all_molecular(registry)
    components = vcat(String[molecular[s]["formula"] for s in neutrals],
                      String[registry.ion[s]["formula"] for s in ions])
    charges = vcat(zeros(Int, length(neutrals)), Int[registry.ion_charge[s] for s in ions])
    ion_lookup = Dict(py_str(row["Species"]) => row for row in registry.ion_epcsaft_table)
    ion_mw = Dict(name => sum(atomic_weights[e] * count for (e, count) in d["elements"]; init=0.0)
                  for (name, d) in registry.ion)
    Mw, segment, sigma, epsilon, n_H, n_e = Float64[], Float64[], Float64[], Float64[], Int[], Int[]
    epsilon_assoc, bondvol = Dict{AssocKey,Float64}(), Dict{AssocKey,Float64}()
    for (index, name) in enumerate(all_species)
        if name in neutrals
            haskey(resolved, name) || throw(KeyError("No resolved molecular PC-SAFT parameters for $name"))
            p = resolved[name]
            default_mw = sum(atomic_weights[e] * c for (e, c) in molecular[name]["elements"]; init=0.0)
        else
            haskey(ion_lookup, name) || throw(KeyError("No resolved ion ePC-SAFT parameters for $name"))
            p = ion_lookup[name]
            default_mw = ion_mw[name]
        end
        Mw_i = _required_float(_first_value(p, ("Mw", "MW", "MW_g_mol", "MW [g/mol]"), default_mw), "Mw", name)
        m_i = _required_float(_first_value(p, ("m", "segment")), "m", name)
        sigma_i = _required_float(_first_value(p, ("sigma", "sigma_A", "sigma [A]")), "sigma", name)
        eps_i = _required_float(_first_value(p, ("epsilon", "u_over_kB_K", "epsilon/kB [K]")), "epsilon", name)
        epsAB = _first_value(p, ("epsilonAB", "epsilon_assoc", "epsilonAB/kB [K]"))
        kapAB = _first_value(p, ("kappaAB", "bondvol"))
        scheme = _first_value(p, ("association_scheme", "Association"), "-")
        push!(Mw, Mw_i)
        push!(segment, m_i)
        push!(sigma, sigma_i)
        push!(epsilon, eps_i)
        if !isna(epsAB) && !isna(kapAB) && isfinite(py_float(epsAB)) && isfinite(py_float(kapAB))
            epsAB, kapAB = py_float(epsAB), py_float(kapAB)
            nH, ne = association_sites(scheme, epsAB, kapAB)
        else
            nH, ne = 0, 0
        end
        push!(n_H, nH)
        push!(n_e, ne)
        component = components[index]
        if nH > 0 && ne > 0
            key = ((component, "H"), (component, "e"))
            epsilon_assoc[key] = epsAB
            bondvol[key] = kapAB
        end
    end
    userlocations = Dict{String,Any}("Mw" => Mw, "segment" => segment, "sigma" => sigma, "epsilon" => epsilon,
                                     "n_H" => n_H, "n_e" => n_e, "epsilon_assoc" => epsilon_assoc,
                                     "bondvol" => bondvol)
    return components, charges, PCSAFT(components; userlocations=userlocations, verbose=false)
end

# ---------------------------------------------------------------- EOS probes
function _balanced_probe_weights(p::ElectrolytePhase)
    weights = Float64[get(p.probe_weights, s, 1.0) for s in p.ions]
    imbalance = _dot(p.z, weights)
    candidates = findall(p.z .* imbalance .< 0)
    if abs(imbalance) > 1e-14
        isempty(candidates) && throw(PhaseError("Both positive and negative ionic carriers are needed for EOS probes"))
        j = candidates[1]
        weights[j] += -imbalance / p.z[j]
    end
    return weights
end

"""
    phi(phase::ElectrolytePhase, n)

Fugacity coefficients for the full component vector `n` (ordered like
`all_species`). The probe composition is made electroneutral (smoothly) and
every component gets a trace amount; the physical amounts used in balances are
unchanged.
"""
function phi(p::ElectrolytePhase, n0)
    if !(n0 isa AbstractVector && length(n0) == length(p.all_species)) || any(<(0), n0)
        throw(PhaseError("Invalid aqueous composition vector"))
    end
    n = collect(float.(n0))
    split, z = length(p.neutrals), p.z
    imbalance = _dot(n[split+1:end], z)
    smoothing = hypot(imbalance, 2e-14)
    pos, neg = findall(>(0), z), findall(<(0), z)
    (isempty(pos) || isempty(neg)) &&
        throw(PhaseError("Aqueous EOS probe requires positive and negative ionic carriers"))
    n[split+pos[1]] += (smoothing - imbalance) / (2 * z[pos[1]])
    n[split+neg[1]] += (smoothing + imbalance) / (-2 * z[neg[1]])
    n[1:split] .+= 1e-14
    n[split+1:end] .+= _balanced_probe_weights(p) .* 1e-14
    phis = fugacity_coefficient(p.eos_model, p.P, p.T, n; phase="liquid")
    if length(phis) != length(n) || any(!isfinite, phis) || any(<=(0), phis)
        throw(EOSFailure("Invalid ePC-SAFT fugacity coefficients"))
    end
    return phis
end

# Intercept of the least-squares line through (x, y) (np.polyfit(x, y, 1)[1]).
function _polyfit_intercept(x, y)
    xm, ym = sum(x) / length(x), sum(y) / length(y)
    slope = sum((x .- xm) .* (y .- ym)) / sum((x .- xm) .^ 2)
    return ym - slope * xm
end

function _infinite_dilution_reference(p::ElectrolytePhase)
    logs = Vector{Vector{Float64}}()
    split = length(p.neutrals)
    scales = (1e-5, 1e-6, 1e-7)
    for scale in scales
        n = zeros(length(p.all_species))
        n[findfirst(==("water"), p.all_species)] = 10.0
        n[split+1:end] .= scale .* _balanced_probe_weights(p)
        push!(logs, log.(phi(p, n)))
    end
    x = sqrt.(collect(scales))
    intercept = [_polyfit_intercept(x, [l[j] for l in logs]) for j in eachindex(p.all_species)]
    ions = split+1:length(p.all_species)
    estimate = logs[3][ions] .- (logs[2][ions] .- logs[3][ions]) .* x[3] ./ (x[2] - x[3])
    maximum(abs.(intercept[ions] .- estimate)) > 0.02 &&
        throw(PhaseError("Ion dilution extrapolation did not converge"))
    return intercept
end

# ---------------------------------------------------------------- G and mu
"""Full component vector (ordered like `all_species`) of the selected amounts `n`."""
function expand(p::ElectrolytePhase, n)
    full = zeros(float(eltype(n)), length(p.all_species))
    full[p.positions] .= n
    return full
end

# mu of the neutrals and ions of the full vector `full` (water first) with fugacity coefficients `phis`.
function _electrolyte_mu(p::ElectrolytePhase, full, phis)
    RT, split, n_ion = R * p.T, length(p.neutrals), length(p.ions)
    neutral, ion = 1:split, length(full)-n_ion+1:length(full)
    q = max.(full, 1e-14)
    x = q ./ sum(q)
    xw = x[1]
    mu_neutral = p.mu0_neutral .+ RT .* log.(x[neutral] .* p.P ./ P0 .* phis[neutral])
    molal = max.(full[ion], 1e-14) ./ (MW_WATER_KG * full[1])
    mu_ions = p.mu0_ion .+ RT .* (log.(molal) .+ log(xw) .+ log.(phis[ion]) .- p.lnphi_infinity[ion])
    G = _dot(full[neutral], mu_neutral) + _dot(full[ion], mu_ions)
    return G, mu_neutral, mu_ions
end

"""G of the full component vector `n` (ordered like `all_species`, water first)."""
function gibbs_full(p::ElectrolytePhase, n)
    if n[1] <= 1e-9
        sum(n[2:end]) > 1e-10 && throw(PhaseError("Aqueous solutes without water"))
        return zero(float(eltype(n)))
    end
    return _electrolyte_mu(p, n, phi(p, n))[1]
end

gibbs(p::ElectrolytePhase, n) = gibbs_full(p, expand(p, n))

# On the optimizer path every probe failure is an EOS failure (cut the step back).
function _phi_for_optimizer(p::ElectrolytePhase, full)
    try
        return phi(p, full)
    catch err
        err isa ValueError && throw(EOSFailure(err.msg))
        rethrow()
    end
end

function gibbs_and_mu(p::ElectrolytePhase, n)
    full = expand(p, n)
    if full[1] <= 1e-9
        sum(full[2:end]) > 1e-10 && throw(PhaseError("Aqueous solutes without water"))
        return 0.0, zeros(length(n))
    end
    G, mu_neutral, mu_ions = _electrolyte_mu(p, full, _phi_for_optimizer(p, full))
    return G, vcat(mu_neutral, mu_ions)[p.positions]
end

mu_at(p::ElectrolytePhase, n) = gibbs_and_mu(p, n)[2]

lnphi(p::ElectrolytePhase, n) = log.(_phi_for_optimizer(p, expand(p, n)))[p.positions]

function hessian(p::ElectrolytePhase, n)
    n[1] <= 1e-9 && return zeros(length(n), length(n))
    # Unselected reference components enter the mole fractions at 1e-14.
    total = sum(max.(n, 1e-14)) + 1e-14 * (length(p.all_species) - length(n))
    return _hessian(p, n, total)
end

"""
    _hessian(phase, n, total)

d mu / d n: analytic ideal part + forward-difference ln(phi) part, symmetrized.
"""
function _hessian(phase, n, total)
    RT = R * phase.T
    m = length(n)
    alive = n .> 1e-14
    ideal = [(i == j && alive[i] ? 1.0 / max(n[i], 1e-14) : 0.0) - (alive[i] & alive[j]) / total
             for i in 1:m, j in 1:m]
    base = lnphi(phase, n)
    jac = zeros(m, m)
    h = 1e-7 * sum(n)
    for j in 1:m
        step = collect(Float64, n)
        step[j] += h
        jac[:, j] = (lnphi(phase, step) .- base) ./ h
    end
    return RT .* (ideal .+ 0.5 .* (jac .+ permutedims(jac)))
end
