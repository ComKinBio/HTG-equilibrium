# Species registry and model loading against the Python app:
#   - tests/reference/data/catalogue.json (as tests/test_reference_catalogue.py)
#   - a user species folder written by the Python app (test/dump_user_species.py)
#   - smoke test of the shipped data/ folder
#   - reading a model from several threads
isdefined(Main, :DataTestSetup) || include("data_setup.jl")

module TestRegistry   # own namespace: the test files share Main

using Main.DataTestSetup
using HTGEquilibrium
using JSON
using HTGEquilibrium: OrderedDict   # OrderedCollections is not a test dependency
using Test

const USER_FOLDER = joinpath(@__DIR__, "data", "user_species")
const USER_EXPECTED = joinpath(@__DIR__, "data", "user_species_python.json")

# Hint texts that were reworded for the case-file workflow (as the Python test).
const REWORDED = Dict("NASA available; no resolved SAFT mapping. Supply species_overrides or EOS data." =>
                      NO_SAFT_REASON)

function mismatches(new, ref)
    bad = [(a, b) for (a, b) in zip(new, ref) if !same(a, b)]
    isempty(bad) || @info "first mismatches" bad[1:min(3, end)]
    return length(bad)
end

@testset "registry" begin

@testset "catalogue.json" begin
    catalogue_ref = reference("catalogue")
    for kind in ("base", "acetone")
        model, ref = MODELS[kind], catalogue_ref[kind]
        @testset "$kind" begin
            new = plain(species(model; include_unavailable=true))
            expected = [merge(r, Dict("Reason" => get(REWORDED, r["Reason"], r["Reason"]))) for r in ref["catalogue"]]
            @test length(new) == length(expected)
            @test mismatches(new, expected) == 0
            for name in ("molecular", "organic", "solid", "ion")
                reg = getfield(model.registry, Symbol(name))
                regref = ref["registries"][name]
                @test collect(keys(reg)) == collect(keys(regref))
                bad = [n for (n, e) in regref
                       if !haskey(reg, n) || !same(Dict(k => get(reg[n], k, nothing) for k in keys(e)), e)]
                isempty(bad) || @info "registry mismatches" name bad[1:min(3, end)]
                @test isempty(bad)
            end
            for table in ("ion_epcsaft_table", "ion_hkf_table")
                @test same(plain(getfield(model.registry, Symbol(table))), ref[table])
            end
            @test same(plain(mineral_inventory(model)), ref["mineral_inventory"])
            # The original default element set; source labels were renamed, so they are not compared.
            columns = ["Species", "Phase", "Elements", "Ready", "Reason"]
            new_pc = [OrderedDict(c => r[c] for c in columns)
                      for r in product_candidates(model, ["Na", "K", "Cs", "Li", "Cl", "Mg"]; temperature_C=400.0)]
            ref_pc = [OrderedDict(c => r[c] for c in columns) for r in ref["ion_product_candidates_400C"]]
            @test same(new_pc, ref_pc)
        end
    end
end

@testset "species() filters and lookups" begin
    model = MODELS["acetone"]
    @test all(r -> r.Ready, species(model))
    @test length(species(model; include_unavailable=true)) == length(model.registry.catalogue)
    @test all(r -> occursin("aqueous", lowercase(r.Phase)), species(model; phase="AQUEOUS", include_unavailable=true))
    hits = species(model; search="METH")
    @test !isempty(hits) && all(r -> occursin("meth", lowercase(r.Species)) ||
                                     occursin("meth", lowercase(something(r.Formula, ""))), hits)
    reg = model.registry
    @test nasa_name(reg, "acetone") == "CH3COCH3" && nasa_name(reg, "carbon") == "CHAR"
    @test (try nasa_name(reg, "nope"); catch e; e.key end) == "No NASA mapping for nope"
    @test HTGEquilibrium.entry(reg, "aqueous", "Na+")["charge"] == 1
    @test (try HTGEquilibrium.entry(reg, "solid", "water"); catch e; e.msg end) |> startswith("solid: 'water' is not registered.")
    @test molecular_weight("CH4") == 12.011 + 4 * 1.008
    @test formula_elements("C2H6O") == OrderedDict("C" => 2.0, "H" => 6.0, "O" => 1.0)
    @test formula_elements("Na+") === nothing
    @test reg.binary_interactions == [(("CO2", "CH4"), 0.065)]
    @test reg.aqueous_reference_ions == ["Na+", "HCO3-", "CO3^2-"]
    @test (try resolve_parameters(reg, model.databases, 600.0, ["ASH"]); catch e; e.key end) ==
          "No PC-SAFT formula mapping for 'ASH'."           # a NASA-only record has no SAFT formula
    @test (try resolve_parameters(reg, model.databases, 600.0, ["butane"]); catch e; e.msg end) ==
          "butane: not found in the molecular database registry"
    @test model.load_args["water_model"] == "SUPCRT92" && model.water_model isa SUPCRT92Water
end

@testset "user species folder (Python fixture)" begin
    expected = JSON.parsefile(USER_EXPECTED; allownan=true)
    model = load_model(FROZEN; user_data=USER_FOLDER)
    reg = model.registry
    T, P = expected["T"], expected["P"]
    @test same(plain(model.user_info), expected["user_info"])
    @test collect(keys(model.user_info)) == collect(keys(expected["user_info"]))
    @test model.user_data == USER_FOLDER
    catalogue = plain(species(model; include_unavailable=true))
    @test length(catalogue) == expected["catalogue_length"]
    touched = Set(r["Species"] for r in expected["catalogue_rows"]) ∪ Set(keys(model.user_info)) ∪ Set(["ethane"])
    @test same([r for r in catalogue if r["Species"] in touched], expected["catalogue_rows"])
    for (name, names) in expected["registry_names"]
        @test collect(keys(getfield(reg, Symbol(name)))) == names
    end
    @test !haskey(reg.molecular, "ethane")       # hidden shipped species
    for (name, e) in expected["entries"]
        @test same(plain(HTGEquilibrium.any_entry(reg, name)), e)
    end
    resolved = resolve_parameters(reg, model.databases, T, collect(keys(expected["saft_parameters"])))
    for (name, e) in expected["saft_parameters"]
        @test same(Dict(k => get(resolved[name], k, nothing) for k in keys(e)), e)
    end
    @test resolved["testanol"]["Source"] == "user_pcsaft.csv" && resolved["methane"]["Database row"] < 0
    @test same(plain(reg.binary_interactions), expected["binary_interactions"])
    @test same(plain(reg.ion_epcsaft_table), expected["ion_epcsaft_table"])
    @test same(plain(reg.ion_hkf_table), expected["ion_hkf_table"])
    names = collect(keys(expected["nasa_g0"]))
    @test allclose_np(nasa_g0(model, names, T), [expected["nasa_g0"][n] for n in names])
    names = collect(keys(expected["solid_g0"]))
    @test allclose_np(solid_g0(model, names, T, P), [expected["solid_g0"][n] for n in names])
    names = collect(keys(expected["hkf_ion_g0"]))
    @test allclose_np(hkf_ion_g0(model, names, 573.15, 240e5), [expected["hkf_ion_g0"][n] for n in names])
    # The shipped model is unchanged by the user layer, and without the layer the original is back.
    plain_model = load_model(FROZEN)
    @test resolve_parameters(plain_model.registry, plain_model.databases, T, ["methane"])["methane"]["m"] ==
          resolved["methane"]["m"] / 1.5
    @test haskey(plain_model.registry.molecular, "ethane") && !haskey(plain_model.registry.solid, "Testite")
    # An empty or missing folder changes nothing.
    empty = load_model(FROZEN; user_data=joinpath(@__DIR__, "data", "no-such-folder"))
    @test plain(empty.registry.catalogue) == plain(plain_model.registry.catalogue)
    @test isempty(empty.user_info)
end

@testset "user folder location" begin
    withenv("HTG_USER_DATA" => "~/somewhere") do
        @test default_user_dir() == joinpath(homedir(), "somewhere")
        @test resolve_user_dir() == joinpath(homedir(), "somewhere")
    end
    withenv("HTG_USER_DATA" => nothing, "XDG_DATA_HOME" => "/xdg") do
        Sys.islinux() && @test default_user_dir() == "/xdg/htg-equilibrium/species"
    end
    @test resolve_user_dir("/a/b") == "/a/b"
    @test resolve_user_dir("/a/b"; enabled=false) === nothing
end

@testset "shipped data/ folder loads" begin
    model = load_model(SHIPPED)
    reg = model.registry
    @test model.data_dir == SHIPPED
    withenv("HTG_DATA_DIR" => SHIPPED) do
        @test default_data_dir() == SHIPPED
    end
    ready = species(model)
    @test !isempty(ready)
    for name in ("water", "hydrogen", "carbon dioxide", "methane", "acetone")
        @test haskey(reg.molecular, name)
    end
    @test Set(reg.aqueous_reference_ions) ⊆ Set(keys(reg.ion))
    molecules = [r.Species for r in ready if r.Phase == "molecular"]
    T = 623.15
    resolved = resolve_parameters(reg, model.databases, T, molecules)
    @test all(d -> all(isfinite, (d["m"], d["sigma"], d["epsilon"], d["Mw"])), values(resolved))
    @test all(isfinite, nasa_g0(model, molecules, T))
    minerals = [n for (n, d) in reg.solid if get(d, "adapter", nothing) == "supcrt_mineral" &&
                                               d["supcrt_mineral"]["Tmax"] >= T]
    @test all(isfinite, solid_g0(model, minerals, T, 240e5))
    @test all(isfinite, hkf_ion_g0(model, collect(keys(reg.ion)), T - 50, 240e5))
end

@testset "a loaded model is read safely from several threads" begin
    model = MODELS["acetone"]
    reg, db = model.registry, model.databases
    molecules = [r.Species for r in species(model) if r.Phase == "molecular"]
    solids = collect(keys(reg.solid))
    temps = collect(range(450.0, 800.0; length=32))
    work(T) = (resolve_parameters(reg, db, T, molecules), nasa_g0(model, molecules, T),
               [HTGEquilibrium.isna(x) ? NaN : x for x in [try solid_g0(model, [s], T, 240e5)[1] catch; NaN end
                                                         for s in solids]])
    serial = [work(T) for T in temps]
    threaded = Vector{Any}(undef, length(temps))
    Threads.@threads for i in eachindex(temps)
        threaded[i] = work(temps[i])
    end
    @test all(isequal(plain(serial[i]), plain(threaded[i])) for i in eachindex(temps))
end

end # registry

end # module TestRegistry
