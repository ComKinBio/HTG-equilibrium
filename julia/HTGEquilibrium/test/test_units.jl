# Display units (Python: tests/test_units.py), plus the Python helpers the
# case files rely on (repr, round, float); with test/data/cases_python.json
# also every format_number / conversion against Python.
isdefined(@__MODULE__, :PY_CASES) || include("python_cases.jl")
using HTGEquilibrium.Units
using HTGEquilibrium: py_round, py_float, py_isclose, py_format_g

const U = HTGEquilibrium.Units

# Temperatures stay well above 0 K: K <-> °C is an offset, so the relative error near 0 K grows.
const UNIT_VALUES = Dict("temperature" => [-200.0, -40.0, 0.0, 0.1, 25.0, 300.0, 374.15, 573.15, 1234.5678],
                         "pressure" => [0.0, 1e-3, 0.1, 1.0, 1.01325, 240.0, 1234.5678, 2.4e7])
close_(a, b, rel=1e-12) = abs(a - b) <= rel * max(abs(a), abs(b)) || a == b

@testset "units (htg/units.py)" begin
    @testset "round trips $unit" for unit in U.UNITS
        for value in UNIT_VALUES[U.quantity(unit)]
            @test close_(U.from_file(U.to_file(value, unit), unit), value)
            @test close_(U.to_file(U.from_file(value, unit), unit), value)
        end
    end

    @testset "conversions" begin
        @test U.to_file(573.15, "K") ≈ 300.0 rtol = 1e-14
        @test U.from_file(300, "K") ≈ 573.15 rtol = 1e-14
        @test U.to_file(24, "MPa") == 240
        @test U.to_file(2.4e7, "Pa") == 240
        @test U.from_file(240, "Pa") == 2.4e7
        @test U.to_file(1, "atm") == 1.01325
        @test U.from_file(1.01325, "atm") == 1
        @test U.to_file(300, "C") == 300 && U.to_file(300, "°C") == 300
        @test U.to_file(240, "bar") == 240
        @test close_(U.pressure_from_bar(U.pressure_to_bar(1, "atm"), "Pa"), 101325)
        @test close_(U.pressure_from_bar(U.pressure_to_bar(1, "MPa"), "atm"), 10 / 1.01325)
        @test isapprox(U.to_file([573.15, 673.15], "K"), [300, 400]; rtol=1e-14)
    end

    @testset "unknown or wrong units" begin
        err = try U.to_file(1, "psi") catch e; e end
        @test err isa U.UnitError && err isa HTGEquilibrium.ValueError
        @test err.msg == "Unknown unit \"psi\". Use one of: °C, K, bar, MPa, Pa, atm."
        @test_throws U.UnitError U.temperature_to_C(1, "bar")
        @test_throws U.UnitError U.pressure_to_bar(1, "K")
    end

    @testset "quantities and labels" begin
        @test [U.quantity(u) for u in U.UNITS] == [fill("temperature", 2); fill("pressure", 4)]
        @test U.label("C") == "°C" && U.label("MPa") == "MPa"
        @test U.field_label("Temperatures", "K") == "Temperatures [K]"
        @test U.field_label("Pressures", "bar") == "Pressures [bar]"
        @test U.field_label("Temperatures", "C") == "Temperatures [°C]"
        @test U.FILE_UNITS == Dict("temperature" => "C", "pressure" => "bar")
    end

    @testset "percent" begin
        @test U.percent_to_fraction(90) == 0.9
        @test U.fraction_to_percent(0.05) == 5
        for percent in [0, 0.5, 3, 33.3333, 90, 100]
            @test close_(U.fraction_to_percent(U.percent_to_fraction(percent)), percent)
        end
    end

    @testset "format_number" begin
        f = U.format_number
        @test f(300.0) == "300"
        @test f(0.1 + 0.2) == "0.3"
        @test f(573.15) == "573.15"
        @test f(299.99999999999994) == "300"
        @test f(1234.5678) == "1234.57"
        @test f(1234567) == "1.23457e6"
        @test f(1.5e-5) == "1.5e-5"
        @test f(-40) == "-40"
        @test f(0) == "0" && f(-0.0) == "0"
        @test f(2.4e7, "Pa") == "2.4e7"
        @test f(101325, "Pa") == "1.01325e5"
        @test f(1e5, "Pa") == "1e5"
        @test f(24, "MPa") == "24"
        @test f(Inf) == "inf" && f(NaN) == "nan"
        @test U.format_list([300, 325.5, 350]) == "300, 325.5, 350"
    end

    @testset "display from file units" begin
        @test U.display(240, "Pa") == "2.4e7"
        @test U.display(240, "MPa") == "24"
        @test U.display(240, "atm") == "236.862"
        @test U.display(300, "K") == "573.15"
        @test U.display(300, "C") == "300"
    end

    @testset "Python helpers" begin
        @test py_repr([1e-5, 1e16, 1.5, 0.1, -0.0, 1e15, 0.0001]) == "[1e-05, 1e+16, 1.5, 0.1, -0.0, 1000000000000000.0, 0.0001]"
        @test py_repr(("a",)) == "('a',)" && py_repr(("a", 1)) == "('a', 1)"
        @test py_repr("it's") == "\"it's\"" && py_repr("a\"b") == "'a\"b'" && py_repr("°C") == "'°C'"
        @test py_repr(OrderedDict("a" => nothing, "b" => true)) == "{'a': None, 'b': True}"
        @test py_float(" 3_50 ") == 350.0 && py_float("inf") == Inf && py_float(true) == 1.0
        @test_throws ArgumentError py_float("0x10")
        @test_throws ArgumentError py_float("1__0")
        @test_throws ArgumentError py_float([1])
        @test py_round(2.675, 2) == 2.67 && py_round(0.125, 2) == 0.12 && py_round(0.375, 2) == 0.38
        @test py_isclose(1.0, 1.0 + 1e-10) && !py_isclose(1.0, 1.0 + 1e-8)
        @test py_isclose(0.0, 1e-11; rel_tol=0, abs_tol=1e-10)
        # Julia's @sprintf("%.6g", 999999.5) is "1000000"; C and Python give "1e+06"
        @test py_format_g(999999.5, 6) == "1e+06" && py_format_g(123456.5, 6) == "123456"
        @test py_format_g(0.0001, 4) == "0.0001" && py_format_g(0.00001, 4) == "1e-05" && py_format_g(91.55, 4) == "91.55"
    end

    if PY_CASES !== nothing
        @testset "against Python: repr, round" begin
            @test all(py_repr(decode(x)) == r for (x, r) in PY_CASES["reprs"])
            @test all(py_round(x, n) === Float64(r) for (x, n, r) in PY_CASES["rounds"])
        end
        @testset "against Python: format_number ($(length(PY_CASES["format_number"])) calls)" begin
            bad = [(x, unit, digits, want, U.format_number(decode(x), unit, digits))
                   for (x, unit, digits, want) in PY_CASES["format_number"]
                   if U.format_number(decode(x), unit, digits) != want]
            isempty(bad) || @error "format_number differs" bad
            @test isempty(bad)
        end
        @testset "against Python: conversions and display" begin
            for (x, unit, to, from, shown) in PY_CASES["conversions"]
                @test_same enc(U.to_file(x, unit)) to "to_file($x, $unit)"
                @test_same enc(U.from_file(x, unit)) from "from_file($x, $unit)"
                @test U.display(x, unit) == shown
            end
        end
    end
end
