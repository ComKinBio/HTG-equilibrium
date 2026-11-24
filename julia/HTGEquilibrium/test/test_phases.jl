# Phase models (phases.jl):
# - tests/reference/data/thermo.json: fluid-phase G and ln phi, electrolyte
#   lnphi_infinity, G and ln phi, with the tolerances of tests/test_reference_thermo.py;
# - the property tests of tests/test_phases.py (mu is the gradient of G, homogeneity);
# - the exact Hessians against central differences of mu, symmetry, Gibbs-Duhem;
# - the water model of the electrolyte phase, errors, threads.
isdefined(Main, :DataTestSetup) || include("data_setup.jl")

module TestPhases   # own namespace: the test files share Main

using Main.DataTestSetup
using HTGEquilibrium
using HTGEquilibrium: ATOMIC_WEIGHTS, phi, expand
using JSON
using Test

kelvin_pascal(state) = (state["T_C"] + 273.15, state["P_bar"] * 1e5)   # as the Python test

# pytest.approx(expected, rel, abs): |actual - expected| <= max(rel |expected|, abs)
approx_py(actual, expected; rel, atol) = abs(actual - expected) <= max(rel * abs(expected), atol)

function fluid_phase(model, label, names, T, P, root; eos_model=nothing)
    reg = model.registry
    if eos_model === nothing
        resolved = resolve_parameters(reg, model.databases, T, names)
        eos_model = build_pcsaft_model(names, resolved, reg.binary_interactions)
    end
    return FluidPhase(label, names, nasa_g0(model, names, T), eos_model, T, P, root)
end

function electrolyte_phase(model, selection, T, P; water_model=SUPCRT92Water())
    reg = model.registry
    neutrals = unique(vcat(reg.aqueous_reference_neutrals, [s for s in selection if !haskey(reg.ion, s)]))
    resolved = resolve_parameters(reg, model.databases, T, neutrals)
    return ElectrolytePhase(water_model, model, selection, T, P, resolved, ATOMIC_WEIGHTS)
end

const MODEL = MODELS["base"]
const STATES = reference("thermo")["states"]

@testset "phases (thermo.json)" begin
    worst = Dict{String,Float64}()
    note!(key, a, b) = for (x, y) in zip(a, b)
        worst[key] = max(get(worst, key, 0.0), reldev(x, y))
    end
    errors = 0
    for state in STATES
        T, P = kelvin_pascal(state)
        @testset "$(state["T_C"]) C, $(state["P_bar"]) bar" begin
            for (label, ref) in state["fluid_phases"]
                names = String.(ref["species"])
                phase = fluid_phase(MODEL, label, names, T, P, ref["root"])
                @test allclose_np(phase.mu0, Float64.(ref["mu0"]))
                for point in ref["points"]
                    n = [Float64(point["amounts"][s]) for s in names]
                    if haskey(point, "error")
                        errors += 1
                        err = try gibbs(phase, n); nothing catch e; e end
                        @test err isa EOSFailure && err isa ValueError
                        # Python's message; x (a list of floats) compared with RTOL.
                        head, x_python = split(chopprefix(point["error"], "ValueError: "), ", x=")
                        @test startswith(err.msg, head * ", x=")
                        x = JSON.parse(split(err.msg, ", x=")[2])
                        @test allclose_np(Float64.(x), Float64.(JSON.parse(x_python)); rtol=RTOL, atol=0.0)
                        continue
                    end
                    G = gibbs(phase, n)
                    @test approx_py(G, point["G"]; rel=1e-10, atol=1e-8)
                    note!("fluid G", G, point["G"])
                    lp = lnphi(phase, n)
                    @test allclose_np(lp, Float64.(point["ln_phi"]); rtol=1e-9, atol=1e-11)
                    note!("fluid ln phi", lp, point["ln_phi"])
                end
            end
            for (label, ref) in get(state, "electrolyte", Dict())
                phase = electrolyte_phase(MODEL, String.(ref["selection"]), T, P)   # the reference used SUPCRT92
                @test phase.all_species == ref["aqueous_species"]
                @test phase.species == ref["selection"]
                @test same(water_dict(phase.water), ref["water"])
                @test allclose_np(phase.mu0_ion, Float64.(ref["mu0_ion"]))
                @test allclose_np(phase.mu0_neutral, Float64.(ref["mu0_neutral"]))
                @test allclose_np(phase.lnphi_infinity, Float64.(ref["lnphi_infinity"]); rtol=1e-9, atol=1e-11)
                note!("electrolyte lnphi_infinity", phase.lnphi_infinity, ref["lnphi_infinity"])
                for point in ref["points"]
                    full = [Float64(get(point["amounts"], s, 0.0)) for s in phase.all_species]
                    if haskey(point, "error")
                        errors += 1
                        @test_throws ValueError gibbs_full(phase, full)
                        continue
                    end
                    G = gibbs_full(phase, full)
                    @test approx_py(G, point["G"]; rel=1e-10, atol=1e-8)
                    note!("electrolyte G", G, point["G"])
                    lp = log.(phi(phase, full))
                    @test allclose_np(lp, Float64.(point["ln_phi"]); rtol=1e-9, atol=1e-11)
                    note!("electrolyte ln phi", lp, point["ln_phi"])
                end
            end
        end
    end
    @test errors == 7   # the reference's EOS failures (all in fluid phases)
    @info "phases vs thermo.json: worst relative deviation" worst
end

# ---------------------------------------------------------------- tests/test_phases.py
const GAS = ["water", "hydrogen", "carbon dioxide", "methane", "nitrogen"]
const AQUEOUS = ["water", "carbon dioxide", "methane", "Na+", "HCO3-"]
const DRY_GAS = ["hydrogen", "carbon dioxide", "methane", "nitrogen"]
const T, P = 300 + 273.15, 240e5

const PHASES = Dict(
    "gas" => (fluid_phase(MODEL, "gas", GAS, T, P, "vapour"), [0.3, 0.05, 0.4, 0.2, 0.05]),
    "aqueous" => (electrolyte_phase(MODEL, AQUEOUS, T, P), [10.0, 0.01, 0.002, 0.003, 0.003]))

@testset "phases: mu is the gradient of G, homogeneity ($name)" for name in ("gas", "aqueous")
    phase, n = PHASES[name]
    G, mu = gibbs_and_mu(phase, n)
    @test approx_py(G, gibbs(phase, n); rel=1e-12, atol=0.0)
    @test G ≈ sum(n .* mu) rtol = 1e-14   # Euler
    for (i, species) in enumerate(phase.species)
        species in ("Na+", "HCO3-") && continue   # a single ion cannot be varied without breaking electroneutrality
        # Step large enough that round-off in G (~1e6 J) stays below the tolerance.
        h = 1e-4 * n[i]
        up, down = copy(n), copy(n)
        up[i] += h
        down[i] -= h
        derivative = (gibbs(phase, up) - gibbs(phase, down)) / (2h)
        @test approx_py(derivative, mu[i]; rel=1e-6, atol=1e-2)
    end
    if name == "aqueous"   # the ions as an electroneutral pair (not in the Python test)
        d = [0.0, 0, 0, 1, 1]
        h = 1e-4 * 0.003
        derivative = (gibbs(phase, n .+ h .* d) - gibbs(phase, n .- h .* d)) / (2h)
        @test approx_py(derivative, sum(mu .* d); rel=1e-6, atol=1e-2)
    end
    @test approx_py(gibbs(phase, 2.5 .* n), 2.5 * gibbs(phase, n); rel=1e-10, atol=0.0)
    @test mu_at(phase, n) == mu
end

# ---------------------------------------------------------------- Hessians
# `hessian` is exact (ForwardDiff Jacobian of mu through Clapeyron; owner, PR #9).
# Checked against central differences of mu, column by column for fluids and along
# electroneutral directions for the electrolyte (the EOS probe makes every composition
# electroneutral, so a single ion cannot be varied alone).
colrel(A, B) = maximum(maximum(abs.(A[:, j] .- B[:, j])) / maximum(abs.(B[:, j])) for j in axes(B, 2))

"""Central differences of mu along the columns of D (steps rel * |amounts moved|)."""
function central_mu(phase, n, D; rel=1e-4)
    return reduce(hcat, begin
                      d = D[:, k]
                      h = rel * minimum(n[i] for i in eachindex(d) if d[i] != 0)
                      (mu_at(phase, n .+ h .* d) .- mu_at(phase, n .- h .* d)) ./ (2h)
                  end for k in axes(D, 2))
end

"""Electroneutral directions: each neutral alone, and the given ion combinations."""
function electroneutral_directions(phase, combinations)
    columns = [Float64.(phase.species .== s) for s in phase.species if !(s in phase.ions)]
    for combination in combinations
        d = zeros(length(phase.species))
        for (s, c) in combination
            d[findfirst(==(s), phase.species)] = c
        end
        @assert sum(d[i] * phase.z[findfirst(==(s), phase.ions)] for (i, s) in enumerate(phase.species)
                    if s in phase.ions) == 0
        push!(columns, d)
    end
    return reduce(hcat, columns)
end

@testset "phases: Hessians" begin
    worst = Dict{String,Float64}()
    note!(key, v) = (worst[key] = max(get(worst, key, 0.0), v))
    state = STATES[2]   # 300 C, 240 bar
    Ts, Ps = kelvin_pascal(state)
    for (label, ref) in state["fluid_phases"]
        names = String.(ref["species"])
        phase = fluid_phase(MODEL, label, names, Ts, Ps, ref["root"])
        for point in ref["points"][1:2]
            n = [max(Float64(point["amounts"][s]), 1e-6) for s in names]
            H = hessian(phase, n)
            @test H == permutedims(H)
            fd = central_mu(phase, n, [Float64(i == j) for i in eachindex(n), j in eachindex(n)])
            @test colrel(H, fd) < 1e-6
            note!("fluid: Hessian vs central differences of mu", colrel(H, fd))
            # Gibbs-Duhem: mu is homogeneous of degree zero, so H n = 0.
            scale = maximum(abs.(H) * n)
            @test maximum(abs.(H * n)) <= 1e-12 * scale
            note!("fluid: |H n| / max |H||n|", maximum(abs.(H * n)) / scale)
        end
    end
    # Electrolyte: the test_phases.py composition and the full_text selection of thermo.json.
    phase, n = PHASES["aqueous"]
    full_text = state["electrolyte"]["full_text"]
    big = electrolyte_phase(MODEL, String.(full_text["selection"]), Ts, Ps)
    nbig = [max(Float64(get(full_text["points"][2]["amounts"], s, 0.0)), 1e-6) for s in big.species]
    # Make nbig electroneutral (as the optimizer's balances do) with Cl- or Na+.
    imbalance = sum(nbig[i] * big.z[findfirst(==(s), big.ions)] for (i, s) in enumerate(big.species) if s in big.ions)
    nbig[findfirst(==(imbalance > 0 ? "Cl-" : "Na+"), big.species)] += abs(imbalance)
    for (ph, amounts, combinations) in (
            (phase, n, [["Na+" => 1, "HCO3-" => 1]]),
            (big, nbig, [["Na+" => 1, "HCO3-" => 1], ["Na+" => 2, "CO3^2-" => 1], ["K+" => 1, "Cl-" => 1],
                         ["Mg2+" => 1, "Cl-" => 2], ["Cs+" => 1, "Li+" => 1, "CO3^2-" => 1]]))
        H = hessian(ph, amounts)
        @test H == permutedims(H)
        Z = electroneutral_directions(ph, combinations)
        fd = central_mu(ph, amounts, Z)
        reduced = Z' * H * Z
        @test colrel(reduced, Z' * fd) < 1e-6
        note!("electrolyte: Z'HZ vs central differences (incl. ion combinations)", colrel(reduced, Z' * fd))
        # H Z itself in the rows of the neutral species (the ion rows carry the probe's
        # split through the symmetrisation; only Z'HZ matters on electroneutral directions).
        # Relative to the column's largest entry: rows of tiny amounts make the central
        # differences noisy in the small entries.
        neutral = findall(s -> !(s in ph.ions), ph.species)
        dev = maximum(maximum(abs.((H * Z)[neutral, k] .- fd[neutral, k])) / maximum(abs.(fd[:, k])) for k in axes(Z, 2))
        @test dev < 1e-6
        note!("electrolyte: (H Z) neutral rows vs central differences", dev)
        # Gibbs-Duhem along the composition itself (an electroneutral direction).
        @test maximum(abs.(Z' * H * amounts)) <= 1e-10 * maximum(abs.(Z') * abs.(H) * amounts)
    end
    @info "Hessians: worst relative deviation" worst
end

# ---------------------------------------------------------------- details and errors
@testset "phases: details" begin
    @test association_sites("2B", NaN, NaN) == (1, 1)
    @test association_sites(" 3b ", NaN, NaN) == (1, 2)
    @test association_sites("-", 2500.0, 0.03) == (1, 1)
    @test association_sites("-", NaN, 0.03) == (0, 0)
    @test association_sites(missing, 2500.0, 0.0) == (0, 0)
    @test EOSFailure <: ValueError && PhaseError <: ValueError

    reg = MODEL.registry
    err = try build_pcsaft_model(["water", "methane"], Dict(), reg.binary_interactions); catch e; e end
    @test err isa KeyError && err.key == "No resolved PC-SAFT parameter set for 'water'."

    gas, n = PHASES["gas"]
    @test gas.name == "gas" && gas.species == GAS && gas.root == "vapour"
    # Empty phase: G = 0 and the pure-component mu (computed once). Pure water has
    # no vapour root at 300 C, 240 bar, so this gas raises EOSFailure there (as in Python).
    message(f) = try f(); "no error" catch e; (typeof(e), e.msg) end
    @test message(() -> gibbs_and_mu(gas, zeros(5)))[1] == EOSFailure
    dry = fluid_phase(MODEL, "gas", DRY_GAS, T, P, "vapour")
    G, mu = gibbs_and_mu(dry, zeros(4))
    @test G == 0.0 && gibbs(dry, zeros(4)) == 0.0
    @test mu == [mu_at(dry, [Float64(i == j) for j in 1:4])[i] for i in 1:4]
    @test gibbs_and_mu(dry, zeros(4))[2] == mu
    @test hessian(dry, zeros(4)) == zeros(4, 4)

    solid = SolidPhase(["Calcite", "Graphite"], [-1.2e6, -2.0e4])
    @test gibbs(solid, [2.0, 3.0]) == 2.0 * -1.2e6 + 3.0 * -2.0e4
    @test gibbs_and_mu(solid, [2.0, 3.0]) == (gibbs(solid, [2.0, 3.0]), [-1.2e6, -2.0e4])

    aq, n = PHASES["aqueous"]
    @test aq.positions == [findfirst(==(s), aq.all_species) for s in AQUEOUS]
    @test expand(aq, n)[aq.positions] == n && sum(expand(aq, n)) == sum(n)
    @test aq.components[1] == "H2O" && aq.charges == vcat(zeros(Int, length(aq.neutrals)),
                                                          [reg.ion_charge[s] for s in aq.ions])
    @test gibbs(aq, zeros(5)) == 0.0 && gibbs_and_mu(aq, zeros(5)) == (0.0, zeros(5))
    @test hessian(aq, zeros(5)) == zeros(5, 5)
    @test message(() -> gibbs(aq, [0.0, 1.0, 0, 0, 0])) == (PhaseError, "Aqueous solutes without water")
    @test message(() -> gibbs_and_mu(aq, [0.0, 1.0, 0, 0, 0])) == (PhaseError, "Aqueous solutes without water")
    @test message(() -> phi(aq, n)) == (PhaseError, "Invalid aqueous composition vector")
    # On the optimizer path every probe failure is an EOSFailure.
    @test message(() -> lnphi(aq, [10.0, -1.0, 0, 0, 0])) == (EOSFailure, "Invalid aqueous composition vector")
    @test message(() -> gibbs_full(aq, vcat(10.0, -1.0, zeros(length(aq.all_species) - 2)))) ==
          (PhaseError, "Invalid aqueous composition vector")
    # Outside the aqueous validity range (450 C, 240 bar: rho < 350 kg/m3).
    w = water_properties(SUPCRT92Water(), 723.15, 240e5)
    @test message(() -> electrolyte_phase(MODEL, AQUEOUS, 723.15, 240e5)) ==
          (PhaseError, "State outside SUPCRT aqueous validity checks: rho=$(py_repr(w.rho)), epsilon=$(py_repr(w.epsilon))")
end

@testset "phases: electrolyte water found by name, not position (owner, PR #9)" begin
    first_water, n = PHASES["aqueous"]
    order = [2, 1, 3, 5, 4]                       # carbon dioxide, water, methane, HCO3-, Na+
    other = electrolyte_phase(MODEL, AQUEOUS[order], T, P)
    @test other.species == AQUEOUS[order] && other.all_species == first_water.all_species
    @test other.water_index == 2 && first_water.water_index == 1
    G, mu = gibbs_and_mu(first_water, n)
    G2, mu2 = gibbs_and_mu(other, n[order])
    @test G2 ≈ G rtol = 1e-14
    @test mu2 ≈ mu[order] rtol = 1e-14
    @test hessian(other, n[order]) ≈ hessian(first_water, n)[order, order] rtol = 1e-12
    # Water present, first selected species (carbon dioxide) absent: not empty.
    m = copy(n[order]); m[1] = 0.0
    @test any(!=(0), hessian(other, m)) && gibbs_and_mu(other, m)[1] != 0
    # No water: empty (zeros), or an error with solutes, whatever the order.
    message(f) = try f(); "no error" catch e; (typeof(e), e.msg) end
    @test hessian(other, [0.0, 0, 0, 0, 0]) == zeros(5, 5) && gibbs_and_mu(other, zeros(5)) == (0.0, zeros(5))
    m = copy(n[order]); m[2] = 0.0
    @test hessian(other, m) == zeros(5, 5)
    @test message(() -> gibbs_and_mu(other, m)) == (PhaseError, "Aqueous solutes without water")
    @test message(() -> gibbs(other, m)) == (PhaseError, "Aqueous solutes without water")
    # Water not selected at all.
    dry = electrolyte_phase(MODEL, ["carbon dioxide", "Na+", "HCO3-"], T, P)
    @test dry.water_index == 0
    @test hessian(dry, [1.0, 0.1, 0.1]) == zeros(3, 3)
    @test message(() -> gibbs_and_mu(dry, [1.0, 0.1, 0.1])) == (PhaseError, "Aqueous solutes without water")
end

@testset "phases: water model of the electrolyte phase" begin
    supcrt, n = PHASES["aqueous"]
    iapws = electrolyte_phase(MODEL, AQUEOUS, T, P; water_model=IAPWS95Water())
    @test supcrt.water_model isa SUPCRT92Water && iapws.water_model isa IAPWS95Water
    # One water model throughout: dielectric constant, validity check and ion g°.
    w = water_properties(IAPWS95Water(), T, P)
    @test iapws.water == w
    @test iapws.mu0_ion == hkf_ion_g0(IAPWS95Water(), MODEL, iapws.ions, T, P)
    @test iapws.mu0_neutral == supcrt.mu0_neutral
    @test iapws.water.epsilon != supcrt.water.epsilon
    @test !isapprox(iapws.mu0_ion, supcrt.mu0_ion; rtol=1e-8)
    @test all(isfinite, iapws.lnphi_infinity) && isfinite(gibbs(iapws, n))
    @test !isapprox(gibbs(iapws, n), gibbs(supcrt, n); rtol=1e-10)
end

# ---------------------------------------------------------------- threads
# The app solves points on threads (M1): each task builds its own phase objects; the
# Clapeyron model of a fluid phase may be shared. Results must equal the serial ones.
function phase_quantities(model, eos_model, state_index)
    state = STATES[state_index]
    Ts, Ps = kelvin_pascal(state)
    ref = state["fluid_phases"]["gas"]
    names = String.(ref["species"])
    gas = fluid_phase(model, "gas", names, Ts, Ps, "vapour"; eos_model=eos_model[state_index])
    out = Any[]
    for point in ref["points"]
        haskey(point, "error") && continue
        n = [Float64(point["amounts"][s]) for s in names]
        push!(out, gibbs_and_mu(gas, n), lnphi(gas, n), hessian(gas, n))
    end
    # Empty phase: the pure-component mu, or the EOSFailure where pure water has no vapour root.
    push!(out, try gibbs_and_mu(gas, zeros(length(names))) catch e; e isa EOSFailure ? e.msg : rethrow() end)
    for (_, eref) in get(state, "electrolyte", Dict())
        aq = electrolyte_phase(model, String.(eref["selection"]), Ts, Ps)
        push!(out, aq.lnphi_infinity, aq.mu0_ion)
        for point in eref["points"]
            n = [Float64(get(point["amounts"], s, 0.0)) for s in aq.species]
            push!(out, gibbs_and_mu(aq, n), lnphi(aq, n), hessian(aq, n))
        end
    end
    return out
end

@testset "phases: threads" begin
    Threads.nthreads() == 1 && @info "Run with JULIA_NUM_THREADS=4 to exercise the threaded phase test on several threads."
    reg = MODEL.registry
    eos_model = map(STATES) do state   # one shared Clapeyron model per state
        names = String.(state["fluid_phases"]["gas"]["species"])
        resolved = resolve_parameters(reg, MODEL.databases, state["T_C"] + 273.15, names)
        build_pcsaft_model(names, resolved, reg.binary_interactions)
    end
    serial = [phase_quantities(MODEL, eos_model, i) for i in eachindex(STATES)]
    jobs = repeat(collect(eachindex(STATES)), 4)
    order = jobs[sortperm(rand(length(jobs)))]
    tasks = [Threads.@spawn (Threads.threadid(), phase_quantities(MODEL, eos_model, i)) for i in order]
    results = fetch.(tasks)
    @test all(isequal(r, serial[i]) for ((_, r), i) in zip(results, order))   # bit for bit
    # One FluidPhase shared by tasks: the pure-component mu cache is filled under a lock.
    shared = fluid_phase(MODEL, "gas", DRY_GAS, T, P, "vapour")
    fresh = fluid_phase(MODEL, "gas", DRY_GAS, T, P, "vapour"; eos_model=shared.eos_model)
    expected = gibbs_and_mu(shared, zeros(4))
    @test all(isequal(r, expected) for r in fetch.([Threads.@spawn(gibbs_and_mu(fresh, zeros(4))) for _ in 1:16]))
    @info "phases: threaded results equal to serial" threads = Threads.nthreads() tasks = length(order) threads_used = length(unique(first.(results)))
end

end # module TestPhases
