# Formula parsing and constants against the Python app (htg/formulas.py).
# Expected values: test/data/formulas_python.json (test/dump_formulas.py).
using HTGEquilibrium
using HTGEquilibrium: ATOMIC_WEIGHTS, elements_text, pysum
using JSON
using HTGEquilibrium: OrderedDict
using Test

const FORMULAS = JSON.parsefile(joinpath(@__DIR__, "data", "formulas_python.json"))

outcome(f, args...) = try
    (ok = f(args...),)
catch e
    e isa FormulaError || rethrow()
    (error = sprint(showerror, e),)
end
pyvalue(x) = x === nothing ? nothing : x   # JSON null -> nothing

@testset "formulas (htg/formulas.py)" begin
    @testset "parse_formula" for (text, want) in FORMULAS["parse"]
        got = outcome(parse_formula, text)
        if haskey(want, "error")
            @test got == (error = want["error"],)
        else
            elements, charge = got.ok
            @test [[k, v] for (k, v) in elements] == want["elements"]
            @test charge == pyvalue(want["charge"])
            @test molar_mass(elements) == want["mass"]
        end
    end
    @testset "formula_summary" for (text, want) in FORMULAS["summary"]
        got = outcome(formula_summary, text)
        @test haskey(want, "error") ? got == (error = want["error"],) : got == (ok = want["ok"],)
    end
    @testset "is_plain_formula" for (text, want) in FORMULAS["plain"]
        @test is_plain_formula(text) == want
    end
    @testset "charge given" for (text, z, split, summary) in FORMULAS["charged"]
        neutral, found = split_charge(text, pyvalue(z))
        @test [neutral, found] == split
        got = outcome(formula_summary, text, pyvalue(z))
        @test haskey(summary, "error") ? got == (error = summary["error"],) : got == (ok = summary["ok"],)
    end
    @testset "hill_formula" for (pairs, want) in FORMULAS["hill"]
        @test hill_formula(OrderedDict(k => Float64(v) for (k, v) in pairs)) == want
    end
    @testset "parse_element_list" for (text, want) in FORMULAS["elements"]
        got = outcome(parse_element_list, text)
        @test haskey(want, "error") ? got == (error = want["error"],) : got == (ok = want["ok"],)
    end
    @testset "pysum (Python 3.12 sum)" begin
        @test pysum([40.078, 31.998, 2.016]) == 74.092          # left fold: 74.09200000000001
        @test pysum(fill(0.1, 10)) == 1.0                       # left fold: 0.9999999999999999
        @test pysum(Float64[]) == 0.0
        @test pysum([1e308, 1e308]) == Inf
        @test pysum([1, 2, 3]) === 6.0
    end
    @test length(ATOMIC_WEIGHTS) == 83
    @test elements_text(OrderedDict("C" => 1.0, "H" => 2.5)) == "C 1 · H 2.5"
end
