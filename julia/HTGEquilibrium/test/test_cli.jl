# The command line (Python: htg/cli.py) against what the Python command line
# prints: test/data/python_cli/ (written by test/dump_cli.py) holds the output,
# stderr and exit code of ~70 Python commands (usage errors, help, errors,
# `htg check` of every case file and of template variants, `htg species`
# variants), the files `htg new` writes, pandas `to_string` of hand-made tables,
# and an `htg run` of three independent points (one failing) with the folder it
# wrote. Text must be equal, except where listed below (help texts and TOML
# syntax errors); the run's tables are compared with compare_cases.jl's rules.
using Test
using HTGEquilibrium
using HTGEquilibrium: OrderedDict, to_string, parse_cli, worker_threads, CASE_TEMPLATE, _pypath
using JSON
using SHA

const CLI_REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CLI_DATA = joinpath(@__DIR__, "data", "python_cli")
const CLI_PY = JSON.parsefile(joinpath(CLI_DATA, "cli_python.json"))
const CLI_USER_SPECIES = joinpath(@__DIR__, "data", "user_species")

module CLICompare   # compare_case: the rules of the example-case comparison
include("compare_cases.jl")
end

"""Run the command line in `dir`: (exit code, stdout, stderr)."""
function run_cli(args; dir=CLI_REPO)
    out, err = IOBuffer(), IOBuffer()
    code = cd(() -> HTGEquilibrium.main(args; out, err), dir)
    return code, String(take!(out)), String(take!(err))
end

"""The usage part of a help text or usage error (up to the first blank line or error line)."""
usage_part(text) = join(Iterators.takewhile(l -> !isempty(l) && !occursin(": error: ", l), split(text, '\n')), '\n')

"""Tokens of a printed table line: numbers within rtol 1e-6 (atol 1e-6), the rest equal."""
function same_tokens(a, b)
    ta, tb = split(a), split(b)
    length(ta) == length(tb) || return false
    for (x, y) in zip(ta, tb)
        x == y && continue
        u, v = tryparse(Float64, x), tryparse(Float64, y)
        (u === nothing || v === nothing) && return false
        abs(u - v) <= 1e-6 + 1e-6 * abs(v) || return false
    end
    return true
end

@testset "command line" begin
    @testset "pandas to_string" begin
        for spec in CLI_PY["tables"]
            columns = spec["columns"]
            value(v) = v == "__nan__" ? NaN : v == "__inf__" ? Inf : v
            if isempty(columns)
                t = DataTable(OrderedDict{String,Any}[]; columns=["Species", "Formula", "Phase"])
            else
                n = length(first(values(columns)))
                t = DataTable([OrderedDict{String,Any}(c => value(v[i]) for (c, v) in columns) for i in 1:n];
                              columns=collect(keys(columns)))
                for (c, dtype) in zip(keys(columns), spec["dtypes"])
                    # pandas keeps None in an object column; DataTable stores missing
                    dtype == "object" && (t.data[c] = Any[v === nothing ? nothing : x
                                                          for (v, x) in zip(columns[c], t[c])])
                    kind = Dict("float64" => :float, "int64" => :int, "bool" => :bool, "str" => :str,
                                "object" => :object)[dtype]
                    @test HTGEquilibrium._pd_kind(t, c) == kind
                end
            end
            @test to_string(t) == spec["text"]
        end
    end

    @testset "argument parsing" begin
        command, v = parse_cli(["run", "a.toml", "b.toml"])
        @test command == "run" && v["case"] == ["a.toml", "b.toml"]
        @test v["workers"] == 1 && v["out"] === nothing && v["no_plots"] == false && v["data"] === nothing
        @test v["user_data"] === nothing && v["no_user_species"] == false
        command, v = parse_cli(["--data", "d", "run", "--work=3", "a.toml", "--no-p", "--out", "o", "--no-u"])
        @test v["workers"] == 3 && v["no_plots"] && v["out"] == "o" && v["data"] == "d" && v["no_user_species"]
        _, v = parse_cli(["run", "--", "-x.toml", "--y"])
        @test v["case"] == ["-x.toml", "--y"]
        _, v = parse_cli(["species", "--phase", "-1", "--search=a=b", "--all"])
        @test v["phase"] == "-1" && v["search"] == "a=b" && v["all"]
        _, v = parse_cli(["gui"])
        @test v["workers"] === nothing && v["port"] === nothing && v["host"] == "127.0.0.1" && v["results"] == "results"
        _, v = parse_cli(["run", "a", "--workers", " 1_0 "])
        @test v["workers"] == 10
        # ambiguous abbreviation (argparse 3.12 message)
        code, out, err = run_cli(["run", "a.toml", "--no"])
        @test code == 2 && endswith(err, "htg run: error: ambiguous option: --no could match --no-plots, --no-user-species\n")
        @test _pypath("results/") == "results" && _pypath("./a//b/./") == "a/b" && _pypath("") == "." &&
              _pypath("/x//y") == "/x/y" && _pypath(".") == "."
    end

    @testset "Python's output: $(join(c["argv"], ' '))" for c in CLI_PY["commands"]
        tmp = mktempdir()
        for (name, text) in CLI_PY["variants"]
            write(joinpath(tmp, "$name.toml"), text)
        end
        args = [replace(a, "{TMP}" => tmp, "{USER_SPECIES}" => CLI_USER_SPECIES) for a in c["argv"]]
        code, out, err = run_cli(args)
        out, err = replace(out, tmp => "{TMP}"), replace(err, tmp => "{TMP}")
        @test code == c["code"]
        help = any(a -> a in ("-h", "--help"), c["argv"]) && c["code"] == 0
        if help   # help texts say what differs in Julia; the usage is Python's
            @test usage_part(out) == usage_part(c["stdout"])
            @test startswith(out, "usage: ")
        elseif haskey(c, "stdout")
            @test out == c["stdout"]
        else
            @test bytes2hex(sha256(out)) == c["stdout_sha256"]
            @test count(==('\n'), out) == c["stdout_lines"]
            @test startswith(out, c["stdout_head"])
        end
        if occursin("bad_syntax", join(c["argv"]))   # TOML syntax errors: Julia's TOML parser's text
            @test startswith(err, "Error: ") && startswith(c["stderr"], "Error: ")
        else
            @test err == c["stderr"]
        end
    end

    @testset "new" begin
        @test read(CASE_TEMPLATE, String) == read(joinpath(CLI_REPO, "src", "htg", "templates", "case_template.toml"), String)
        dir = mktempdir()
        mkdir(joinpath(dir, "sub.dir"))
        for py in CLI_PY["new"]
            code, out, err = run_cli(py["argv"]; dir)
            @test (code, out, err) == (py["code"], py["stdout"], py["stderr"])
        end
        @test read(joinpath(dir, "my_new_case.toml"), String) == CLI_PY["new"][1]["file"]
        @test read(joinpath(dir, "sub.dir", "other"), String) == CLI_PY["new"][2]["file"]
        @test load_case(joinpath(dir, "my_new_case.toml")).name == "my_new_case"
        code, out, err = run_cli(["new", "my_new_case.toml"]; dir)
        @test code == 1 && err == "my_new_case.toml already exists; choose another name.\n"
    end

    @testset "plot and gui are not available yet" begin
        for args in (["plot", "results/x"], ["gui"], ["gui", "--browser", "--port", "8080"])
            code, out, err = run_cli(args)
            @test code == 2 && isempty(out)
            @test startswith(err, "Error: htg $(args[1]) is not available in the Julia version yet")
        end
    end

    @testset "worker threads" begin
        io = IOBuffer()
        n = Threads.nthreads(:default)
        @test worker_threads(n, io) == n && isempty(take!(io))
        @test worker_threads(n + 3, io) == n
        @test occursin("Julia has $n thread(s); using $n worker(s) instead of $(n + 3)", String(take!(io)))
    end

    @testset "run (three independent points, one failing)" begin
        py = CLI_PY["run"]
        dir = mktempdir()
        mkpath(joinpath(dir, "comparison"))
        cp(joinpath(CLI_REPO, "cases", "comparison", "yakaboylu2013_fig7_8.csv"),
           joinpath(dir, "comparison", "yakaboylu2013_fig7_8.csv"))
        write(joinpath(dir, "cli_run_small.toml"), py["case"])
        code, out, err = run_cli(py["argv"]; dir)
        @test code == py["code"] == 1
        @test isempty(err)
        julia_lines = [l for l in split(out, '\n') if !startswith(l, "Finished in ")]
        python_lines = [l for l in split(py["stdout"], '\n') if !startswith(l, "Finished in ")]
        @test length(julia_lines) == length(python_lines)
        table_start = findfirst(==("cli_run_small:"), python_lines)
        for (i, (j, p)) in enumerate(zip(julia_lines, python_lines))   # exact before the summary, then tolerances
            ok = i < table_start ? j == p : same_tokens(j, p)
            ok || @info "run output differs" julia = j python = p
            @test ok
        end
        rep = CLICompare.compare_case("cli_run_small", joinpath(dir, "out"), joinpath(CLI_DATA, "run"))
        isempty(rep.problems) || @info "run results differ" rep.problems
        @test isempty(rep.problems) && rep.points == 4 && rep.matched == 4

        # the same on two worker threads, with the plot note
        if Threads.nthreads(:default) >= 2 && get(ENV, "HTG_FULL_TESTS", "") == "1"
            args = ["run", "cli_run_small.toml", "--out", "out2", "--workers", "2", "--no-user-species"]
            code, out, err = run_cli(args; dir)
            @test code == 1 && isempty(err)
            @test count("Plots are not written by the Julia version yet", out) == 1
            @test occursin("Starting 2 worker threads", out)
            @test occursin("cli_run_small: 300 °C, 240 bar: solved (start: feasible LP)", out)
            rep = CLICompare.compare_case("cli_run_small", joinpath(dir, "out2"), joinpath(CLI_DATA, "run"))
            @test isempty(rep.problems) && rep.matched == 4
        end
    end
end
