# Standard chemical potentials and SAFT parameters against
# tests/reference/data/thermo.json, with the tolerances of
# tests/test_reference_thermo.py. The phase parts of that file (Gibbs
# energies, ln phi, lnphi_infinity) belong to the phases port.
isdefined(Main, :DataTestSetup) || include("data_setup.jl")

module TestStandardStates   # own namespace: the test files share Main

using Main.DataTestSetup
using HTGEquilibrium
using Test

const PARAM_KEYS = ("Formula", "Source", "Database row", "Model", "RefID", "m", "sigma", "epsilon",
                    "epsilonAB", "kappaAB", "association_scheme", "Mw")

kelvin_pascal(state) = (state["T_C"] + 273.15, state["P_bar"] * 1e5)   # as the Python test

@testset "standard states (thermo.json)" begin
    states = reference("thermo")["states"]
    model = MODELS["base"]
    reg, nasa = model.registry, model.databases.nasa
    worst = Dict{String,Float64}()
    note!(key, a, b) = (a isa Real && b isa Real) && (worst[key] = max(get(worst, key, 0.0), reldev(a, b)))
    for state in states
        T, P = kelvin_pascal(state)
        @testset "$(state["T_C"]) C, $(state["P_bar"]) bar" begin
            @test sort(collect(keys(all_molecular(reg)))) == collect(keys(state["molecular_g0"]))
            bad = String[]
            for (name, expected) in state["molecular_g0"]
                actual = outcome(() -> nasa7_properties(nasa, nasa_name(reg, name), T)[4])
                same(actual, expected) || push!(bad, name)
                note!("molecular_g0", actual, expected)
            end
            @test isempty(bad)

            @test sort(collect(keys(reg.solid))) == collect(keys(state["solid_g0"]))
            bad = String[]
            for (name, expected) in state["solid_g0"]
                actual = outcome(() -> solid_g0(reg, nasa, [name], T, P)[1])
                same(actual, expected) || push!(bad, "$name: $actual vs $expected")
                note!("solid_g0", actual, expected)
            end
            isempty(bad) || @info "solid_g0 mismatches" bad
            @test isempty(bad)

            bad = String[]
            for (name, expected) in state["saft_parameters"]
                actual = outcome() do
                    d = resolve_parameters(reg, model.databases, T, [name])[name]
                    Dict(k => get(d, k, nothing) for k in PARAM_KEYS)
                end
                same(actual, expected) || push!(bad, name)
                if actual isa AbstractDict
                    for k in ("m", "sigma", "epsilon", "Mw")
                        note!("saft_parameters", actual[k], expected[k])
                    end
                end
            end
            isempty(bad) || @info "saft_parameters mismatches" bad
            @test isempty(bad)

            # Standard states of the phases (the phase models themselves are a later step).
            for (label, ref) in state["fluid_phases"]
                mu0 = nasa_g0(reg, nasa, ref["species"], T)
                @test allclose_np(mu0, Float64.(ref["mu0"]))
                foreach(((a, b),) -> note!("fluid mu0", a, b), zip(mu0, ref["mu0"]))
                @test resolve_parameters(reg, model.databases, T, ref["species"]) isa AbstractDict
            end
            for (label, ref) in get(state, "electrolyte", Dict())
                selection = ref["selection"]
                # ElectrolytePhase: reference neutrals and ions first, then the selected species.
                neutrals = unique(vcat(reg.aqueous_reference_neutrals, [s for s in selection if !haskey(reg.ion, s)]))
                ions = unique(vcat(reg.aqueous_reference_ions, [s for s in selection if haskey(reg.ion, s)]))
                @test vcat(neutrals, ions) == ref["aqueous_species"]
                w = water_properties(model, T, P)
                @test same(Dict(k => getfield(w, Symbol(k)) for k in keys(ref["water"])), ref["water"])
                mu0_ion = hkf_ion_g0(model, ions, T, P)
                @test allclose_np(mu0_ion, Float64.(ref["mu0_ion"]))
                foreach(((a, b),) -> note!("mu0_ion", a, b), zip(mu0_ion, ref["mu0_ion"]))
                mu0_neutral = nasa_g0(model, neutrals, T)
                @test allclose_np(mu0_neutral, Float64.(ref["mu0_neutral"]))
                foreach(((a, b),) -> note!("mu0_neutral", a, b), zip(mu0_neutral, ref["mu0_neutral"]))
                @test resolve_parameters(reg, model.databases, T, neutrals) isa AbstractDict
            end
        end
    end
    @info "thermo.json: worst relative deviation" worst
end

@testset "standard states: details" begin
    model = MODELS["base"]
    reg, nasa = model.registry, model.databases.nasa
    calcite = reg.solid["Calcite"]["supcrt_mineral"]
    p = supcrt_mineral_properties(calcite, "Calcite", 298.15, 1e5)
    @test p.G ≈ (calcite["Hf"] - 298.15 * calcite["Sr"]) * 4.184 && p.V_m3_per_mol == calcite["Vr"] * 1e-6
    @test (try supcrt_mineral_properties(calcite, "Calcite", 200.0, 1e5); catch e; e.msg end) ==
          "Calcite: T=200 K outside supported [298.15, $(HTGEquilibrium._g(calcite["Tmax"]))] K"
    @test (try supcrt_mineral_properties(calcite, "Calcite", 500.0, 0.0); catch e; e.msg end) ==
          "Invalid mineral pressure"
    @test (try nasa_g0(reg, Dict{String,Any}(), ["water"], 500.0); catch e; e.key end) ==
          "NASA species 'H2O' for 'water' is not present in nasa_db."
    # The registry's HKF records are the ones test_water.jl uses (test/data/hkf_ions.toml).
    fixture = HTGEquilibrium.TOML.parsefile(joinpath(@__DIR__, "data", "hkf_ions.toml"))
    @test hkf_records(reg.ion_hkf_table, collect(keys(reg.ion))) ==
          [HKFRecord(fixture[s]) for s in keys(reg.ion)]
    @test water_sigma(298.15) ≈ 2.7927 + 10.11 * exp(-0.01775 * 298.15) - 1.417 * exp(-0.01146 * 298.15)
end

end # module TestStandardStates
