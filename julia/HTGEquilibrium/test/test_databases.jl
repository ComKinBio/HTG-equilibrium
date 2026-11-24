# Database readers (databases.jl) and the Python-compatibility helpers they use
# (pycompat.jl). The full content of the readers is checked through the golden
# data in test_registry.jl and test_standard_states.jl.
isdefined(Main, :DataTestSetup) || include("data_setup.jl")

module TestDatabases   # own namespace: the test files share Main

using Main.DataTestSetup
using HTGEquilibrium
using HTGEquilibrium: py_float, py_int, py_repr, py_splitlines, parse_toml, cell, nrows, findtext,
                      findpath, finddescendant, iterxml, getattribute, nasa_float, casefold
using HTGEquilibrium: OrderedDict   # OrderedCollections is not a test dependency
using TOML
using Test

@testset "databases" begin

@testset "Python helpers" begin
    @test py_float(" 1_000.5 ") == 1000.5
    @test py_float("1e3") == 1000.0 && py_float("-inf") == -Inf && isnan(py_float("NaN"))
    @test_throws ArgumentError py_float("1,5")
    @test (try py_float(""); catch e; e.msg end) == "could not convert string to float: ''"
    @test (try py_float(nothing); catch e; e.msg end) ==
          "float() argument must be a string or a real number, not 'NoneType'"
    @test py_int("12") == 12 && py_int(3.9) == 3
    @test (try py_int("1.0"); catch e; e.msg end) == "invalid literal for int() with base 10: '1.0'"
    @test [py_repr(x) for x in (1e-5, 1e16, 100.0, 0.1, -2.5e-7, 123456.789, 1.0e15)] ==
          ["1e-05", "1e+16", "100.0", "0.1", "-2.5e-07", "123456.789", "1000000000000000.0"]
    @test py_repr("it's") == "\"it's\"" && py_repr("a") == "'a'"
    @test py_splitlines("a\r\nb\rc\n\nd") == ["a", "b", "c", "", "d"]
    @test py_splitlines("a\r\nb\n"; keepends=true) == ["a\r\n", "b\n"]
    @test casefold("Straße CO₂") == "strasse co₂"
    @test normalize_name(" CO₂ ") == "co2" && normalize_name(missing) == "" && normalize_name(NaN) == ""
    @test normalize_ion_name("CO₃ ²⁻") == "co3²-" && normalize_ion_name(" Na ⁺") == "na+"
    @test nasa_float(" 1.5D+02") == 150.0 && isnan(nasa_float("   "))
end

@testset "TOML in file order (as tomllib)" begin
    text = """
    # comment
    top = 1
    [b]
    z = { y = 1, "x x" = 2, w = { v = 3, u = 4 } }
    a = [1, [2, 3], "]"]
    [a."q.q"]
    k2 = '''multi
    line'''
    k1 = "s#"
    [[arr]]
    n = 2
    m = 1
    [[arr]]
    m = 3
    n = 4
    """
    d = parse_toml(text)
    @test d == TOML.parse(text)
    @test collect(keys(d)) == ["top", "b", "a", "arr"]
    @test collect(keys(d["b"])) == ["z", "a"]
    @test collect(keys(d["b"]["z"])) == ["y", "x x", "w"]
    @test collect(keys(d["b"]["z"]["w"])) == ["v", "u"]
    @test collect(keys(d["a"]["q.q"])) == ["k2", "k1"]
    @test [collect(keys(x)) for x in d["arr"]] == [["n", "m"], ["m", "n"]]   # per table, as tomllib
    config = read_toml(joinpath(FROZEN, "species.toml"))
    @test collect(keys(config)) == ["nasa_names", "molecular", "organic", "solid", "ion", "aqueous_reference",
                                    "binary_interaction"]
    @test collect(keys(config["molecular"]))[1:3] == ["water", "hydrogen", "carbon dioxide"]
    @test collect(keys(config["ion"])) == ["Na+", "HCO3-", "CO3^2-"]
    @test collect(keys(config["solid"]["Lignin1"]["elements"])) == ["C", "H", "O"]
    @test collect(keys(config["aqueous_reference"]["probe_weights"])) == ["Na+", "HCO3-", "CO3^2-"]
end

db = MODELS["base"].databases

@testset "SAFT CSV tables as pandas reads them" begin
    @test nrows(db.pcsaft) == 115 && nrows(db.epcsaft) == 222
    @test db.pcsaft.columns[1] == "record_type" && db.epcsaft.columns[1] == "record_type"   # BOM dropped
    @test length(db.pcsaft.columns) == 33 && length(db.epcsaft.columns) == 23
    @test db.pcsaft.index == 0:114
    # pandas dtypes (pandas 3.1.1): numeric columns
    @test db.epcsaft.numeric == Set(["Charge", "Year", "lij"])
    @test db.pcsaft.numeric == Set(["Charge", "Dipole_D", "Dipole_fraction", "Ionization_energy_eV", "K1", "K2",
                                    "MW_g_mol", "Quadrupole_DA", "Quadrupole_fraction", "Source_extra_value",
                                    "Year", "epsilonAB_over_kB_K", "kappaAB", "kij", "m", "sigma_A", "u_over_kB_K"])
    @test cell(db.epcsaft, "Value", 2) == "T-dependent: 2.7927 + 10.11 exp(-0.01775T) - 1.417 exp(-0.01146T)"
    @test cell(db.epcsaft, "Value", 1) == "1.2047"                    # text column: kept as text
    @test cell(db.pcsaft, "m", 2) === 1.6069
    @test isnan(cell(db.pcsaft, "epsilonAB_over_kB_K", 2))
    @test ismissing(cell(db.epcsaft, "Species2", 1))
    @test db.pc_columns == OrderedDict("m" => "m", "sigma" => "sigma_A", "epsilon" => "u_over_kB_K",
                                       "epsilonAB" => "epsilonAB_over_kB_K", "kappaAB" => "kappaAB",
                                       "association" => "Association_scheme")
end

@testset "supcrt92.xml as ElementTree reads it" begin
    species = iterxml(db.supcrt, "Species")
    @test length(species) == 1553
    @test db.supcrt.tag == "Database"
    first_species = species[1]
    @test findtext(first_species, "Name") == "Akermanite"
    @test findtext(first_species, "Missing") === nothing
    hkf = findpath(first_species, "Thermo/HKF")
    @test findtext(hkf, "Gf") == "-879362.0" && getattribute(findpath(hkf, "Vr"), "units") == "cm3/mol"
    @test finddescendant(first_species, "Date").text == "15.Mar.90"
    @test count(n -> findtext(n, "Type") == "Aqueous", species) == 1338
end

@testset "NASA-7 parser" begin
    @test length(db.nasa) == 229                      # as Python read_nasa_database
    @test count(r -> r.phase == "S", values(db.nasa)) == 41 && count(r -> r.phase == "", values(db.nasa)) == 187
    ch4 = db.nasa["CH4"]
    @test ch4.phase == "" && ch4.elements == OrderedDict("C" => 1.0, "H" => 4.0)
    @test (ch4.T_low, ch4.T_mid, ch4.T_high) == (300.0, 1000.0, 5000.0)
    @test length(ch4.low) == 7 && length(ch4.high) == 7
    @test parse_nasa_header_elements(rpad("X", 24) * "C   1H   4O   1N   0G") == OrderedDict("C" => 1.0, "H" => 4.0, "O" => 1.0)
    @test parse_nasa_header_elements(rpad("X", 24) * "CL  2NA  1") == OrderedDict("Cl" => 2.0, "Na" => 1.0)
    # A record of the user fixture, with \r\n line ends and lines around it.
    text = read(joinpath(@__DIR__, "data", "user_species", "user_thermo.txt"), String)
    lines = py_splitlines(replace("junk\n" * text, "\n" => "\r\n"); keepends=true)
    records = read_nasa_lines(lines)
    @test collect(keys(records)) == ["testanol", "testol", "mychar"]
    @test records["testanol"].low == db.nasa["CH4"].low && records["mychar"].phase == "S"
    @test isempty(read_nasa_lines(["THERMO\n", "END\n"]))
    bad = [rpad("BAD", 45) * "G" * rpad("  300.x", 33) * "1\n", rpad("", 79) * "2\n", rpad("", 79) * "3\n",
           rpad("", 79) * "4\n"]
    @test (try read_nasa_lines(bad); catch e; e.msg end) ==
          "Could not parse NASA temperature limits for species 'BAD'.\nHeader:\n" * rstrip(bad[1], '\n')
    # nasa7_properties: errors as Python
    @test (try nasa7_properties(db.nasa, "CH4", 100.0); catch e; e.msg end) ==
          "Temperature 100.00 K is outside the NASA range for CH4: 300.00–5000.00 K."
    @test (try nasa7_properties(db.nasa, "nope", 300.0); catch e; e.key end) == "NASA species not found: nope"
    @test_throws ArgumentError nasa7_properties(db.nasa, "CH4", 0.0)
    # Python: (41.19355675907371, -70999.48805231652, 197.25406447566482, -149901.11384258245)
    @test all(isapprox.(nasa7_properties(db.nasa, "CH4", 400.0),
                        (41.19355675907371, -70999.48805231652, 197.25406447566482, -149901.11384258245);
                        rtol=1e-14))
end

@testset "load_databases" begin
    @test (try load_databases(joinpath(@__DIR__, "nowhere")); catch e; e.msg end) ==
          "nasa database not found: $(joinpath(@__DIR__, "nowhere", "Thermdat_WasteMech.txt"))"
    paths = Dict(k => joinpath(FROZEN, v) for (k, v) in HTGEquilibrium.DATABASE_FILES)
    db2 = load_databases(; paths)
    @test collect(keys(db2.nasa)) == collect(keys(db.nasa))
    @test db2.paths == OrderedDict(k => paths[k] for k in ("nasa", "pcsaft", "epcsaft", "supcrt"))
end

end # databases

end # module TestDatabases
