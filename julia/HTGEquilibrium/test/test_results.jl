# Result tables, CSV files and chosen results (Python: htg/results.py, and the
# saved-run loader htg/gui/runs.py `load_saved`).
#
# test/data/python_results/ (written by test/dump_results.py) holds a run saved
# by the Python app from hand-made Results (spec.json) and what Python reads
# back from it (python_read.json). Julia builds the same Results and must write
# the same files; it must read Python's files as pandas does; and files Julia
# writes must read back to the same tables.
using Test
using HTGEquilibrium
using HTGEquilibrium: OrderedDict, Pair2, save
using JSON

const PY_RESULTS = joinpath(@__DIR__, "data", "python_results")
const PY_RUN = joinpath(PY_RESULTS, "run")
const RESULT_FILES = ("summary.csv", "species.csv", "phases.csv", "alternatives.csv")
const RESULTS_TEMPLATE = joinpath(@__DIR__, "..", "..", "..", "src", "htg", "templates", "case_template.toml")

_pairs(triples) = OrderedDict{Pair2,Float64}((String(p), String(s)) => Float64(v) for (p, s, v) in triples)
_dict(d) = OrderedDict{String,Any}(String(k) => v for (k, v) in d)

"""The Result of spec.json's point (as dump_results.py builds it)."""
function spec_result(p)
    amounts = _pairs(p["amounts"])
    species_data = OrderedDict{Pair2,OrderedDict{String,Any}}(
        (String(ph), String(s)) => OrderedDict{String,Any}("formula" => d["formula"], "elements" => _dict(d["elements"]),
                                                           "charge" => d["charge"])
        for (ph, s, d) in p["species_data"])
    alternatives = [Alternative(a["label"], a["gibbs"], _pairs(a["amounts"]), a["n_starts"], a["start_labels"];
                                phase_gibbs=_dict(a["phase_gibbs"])) for a in p["alternatives"]]
    return Result(temperature_C=p["T_C"], pressure_bar=p["P_bar"], amounts=amounts, initial_amounts=copy(amounts),
                  gibbs_J_per_mol_C=p["gibbs_J_per_mol_C"],
                  phase_gibbs_J_per_mol_C=OrderedDict{String,Float64}(_dict(p["phase_gibbs_J_per_mol_C"])),
                  balance_residual=p["balance_residual"], charge_residual=p["charge_residual"],
                  iterations=p["iterations"], eos_trial_errors=String[], feed_C_mol_per_g_dry=p["feed_C_mol_per_g_dry"],
                  elements=collect(String, p["elements"]), feed_atoms=Float64.(p["feed_atoms"]),
                  product_atoms=Float64.(p["feed_atoms"]), gas_carbon_fraction=p["gas_carbon_fraction"],
                  carbon_bound_violation=p["carbon_bound_violation"], phase_rules=_dict(p["phase_rules"]),
                  optimizer=_dict(p["optimizer"]), atomic_weights=Dict{String,Float64}(_dict(p["atomic_weights"])),
                  species_data=species_data, alternatives=alternatives)
end

function spec_sweep()
    spec = JSON.parsefile(joinpath(PY_RESULTS, "spec.json"))
    case = replace(load_case(RESULTS_TEMPLATE); temperature_C=spec["case_temperature_C"], pressure_bar=240.0,
                   aqueous_cutoff_C=nothing)
    points = [(p["T_C"], p["P_bar"], haskey(p, "error") ? String(p["error"]) : spec_result(p)) for p in spec["points"]]
    return SweepResult(case, points)
end

"""pandas' dtype name of a Julia table column ("float64", "int64", "bool", "str", "object")."""
function pandas_dtype(t, c)
    values = t[c]
    c in t.numeric && return all(v -> v isa Int, values) ? "int64" : "float64"
    !isempty(values) && all(v -> v isa Bool, values) && return "bool"
    all(v -> v === missing || v isa String, values) && return "str"
    return "object"
end

jsonval(v) = v === missing || (v isa AbstractFloat && isnan(v)) ? nothing : v
"""A Julia table equals Python's (python_read.json): columns, dtypes, values (floats rtol 1e-12)."""
function same_as_python(t, py)
    t.columns == py["columns"] || return false
    [pandas_dtype(t, c) for c in t.columns] == py["dtypes"] || return false
    nrows(t) == length(py["rows"]) || return false
    for (i, row) in enumerate(py["rows"]), (j, c) in enumerate(t.columns)
        a, b = jsonval(t[c][i]), row[j]
        ok = a isa AbstractFloat && b isa Real ? isapprox(a, b; rtol=1e-12) : isequal(a, b)
        ok || (@info "differs" c i a b; return false)
    end
    return true
end

"""Cells of two CSV texts equal (floats within rtol 1e-14 where the text differs)."""
function same_csv(a::String, b::String)
    a == b && return true
    ra, rb = split(a, '\n'), split(b, '\n')
    length(ra) == length(rb) || return false
    for (x, y) in zip(ra, rb)
        x == y && continue
        fx, fy = split(x, ','), split(y, ',')
        length(fx) == length(fy) || return false
        for (u, v) in zip(fx, fy)
            u == v && continue
            pu, pv = tryparse(Float64, u), tryparse(Float64, v)
            (pu !== nothing && pv !== nothing && occursin(r"[.e]", u) && isapprox(pu, pv; rtol=1e-14)) ||
                (@info "CSV cell differs" u v; return false)
        end
    end
    return true
end

"""Values equal after a CSV round trip: an empty string reads back as missing (as in pandas)."""
roundtrip_equal(x, y) = isequal(x, y) || ((x === missing || x == "" || (x isa Float64 && isnan(x))) &&
                                          (y === missing || y == "" || (y isa Float64 && isnan(y))))

@testset "results" begin
    @testset "tables written as pandas writes them" begin
        t = DataTable([OrderedDict{String,Any}("T_C" => 300.0, "P_bar" => 240.0, "status" => "failed",
                                               "error" => "x, \"y\"\nz"),
                       OrderedDict{String,Any}("T_C" => 350.0, "P_bar" => 240.0, "status" => "solved",
                                               "gibbs" => -1.5e-7, "it" => 23, "aq" => false, "aq2" => nothing,
                                               "d" => NaN, "n" => 1, "s" => "ipopt", "notes" => "", "error" => "")])
        # the same table in pandas 3.1.1 (pd.DataFrame(rows).to_csv(index=False))
        @test csv_text(t) == "T_C,P_bar,status,error,gibbs,it,aq,aq2,d,n,s,notes\n" *
                             "300.0,240.0,failed,\"x, \"\"y\"\"\nz\",,,,,,,,\n" *
                             "350.0,240.0,solved,,-1.5e-07,23.0,False,,,1.0,ipopt,\n"
        @test csv_text(DataTable()) == "\n"
        @test csv_text(DataTable(OrderedDict{String,Any}[]; columns=["a", "b"])) == "a,b\n"
        @test csv_text(DataTable([OrderedDict{String,Any}("a" => 1e16, "b" => 1e-5, "c" => 0.1 + 0.2, "d" => -0.0,
                                                          "e" => 123456789012345678.0, "f" => true)])) ==
              "a,b,c,d,e,f\n1e+16,1e-05,0.30000000000000004,-0.0,1.2345678901234568e+17,True\n"
        @test csv_text(DataTable([OrderedDict{String,Any}("a" => nothing), OrderedDict{String,Any}("a" => nothing)])) ==
              "a\n\"\"\n\"\"\n"
        @test csv_text(DataTable([OrderedDict{String,Any}("a" => " x"), OrderedDict{String,Any}("a" => "")])) ==
              "a\n x\n\"\"\n"
        @test csv_text(DataTable([OrderedDict{String,Any}("a" => 1.0, "b" => 3)])) == "a,b\n1.0,3\n"
        # concat re-infers: an int column next to a float column becomes float
        c = concat([DataTable([OrderedDict{String,Any}("x" => 1)]), DataTable([OrderedDict{String,Any}("x" => 2.5)])])
        @test c["x"] == [1.0, 2.5] && csv_text(c) == "x\n1.0\n2.5\n"
    end

    sweep = spec_sweep()
    @testset "same files as the Python app" begin
        dir = mktempdir()
        save(sweep, dir)
        @test sort(readdir(dir)) == sort([f for f in readdir(PY_RUN) if f != "choices.toml"])
        for f in RESULT_FILES
            @test same_csv(read(joinpath(dir, f), String), read(joinpath(PY_RUN, f), String))
        end
        @test read(joinpath(dir, "case.toml"), String) == read(joinpath(PY_RUN, "case.toml"), String)
        # choices.toml: the same text; a choice of A is left out
        @test write_choices(dir, Dict((400.0, 240.0) => "B", (450, 240) => "A", (300.0, 240.0) => "Z")) ==
              joinpath(dir, CHOICES_FILE)
        @test read(joinpath(dir, CHOICES_FILE), String) == read(joinpath(PY_RUN, CHOICES_FILE), String)
        @test write_choices(dir, Dict((400.0, 240.0) => "A")) === nothing && !isfile(joinpath(dir, CHOICES_FILE))
        @test read_choices(dir) == OrderedDict{Tuple{Float64,Float64},String}()
    end

    @testset "tables of one result and of the sweep" begin
        r = results(sweep)[1]
        s = summary(r)
        @test collect(keys(s)) == ["T_C", "P_bar", "status", "gibbs_J_per_mol_C", "gas_carbon_percent",
                                   "balance_residual", "charge_residual", "iterations", "aqueous_removed", "solver",
                                   "descent_RT", "n_results", "carbon_bound_violation", "notes", "error"]
        @test s["status"] == "solved" && s["aqueous_removed"] == false && s["n_results"] == 1
        t = species_table(r)
        @test nrows(t) == length(r.amounts) && t["group"][7] == "aqueous ions"
        @test !isnan(t["molality_mol_per_kg_water"][7]) && isnan(t["molality_mol_per_kg_water"][1])
        p = phase_summary(r)
        @test p["phase"] == ["gas", "aqueous", "organic_liquid", "solid"]
        @test isapprox(sum(p["carbon_percent"]), 100 * sum(v * get(r.species_data[k]["elements"], "C", 0)
                                                           for (k, v) in r.amounts) / r.feed_atoms[1])
        @test length(results(sweep)) == 3
        st = summary(sweep)
        @test st["status"] == ["solved", "failed", "solved", "solved"]
        @test "iterations" in st.numeric && st["iterations"][2] isa Float64 && isnan(st["iterations"][2])
        @test st["notes"][3] == "aqueous carbon bound dropped"
        alt = alternatives_table(sweep)
        @test alt.columns[1:8] == ["T_C", "P_bar", "result", "gibbs_J_per_mol_C", "n_starts", "phase", "species",
                                   "mol_per_mol_C"]
        @test unique(alt["result"]) == ["A", "B", "C"] && nrows(alt) == 5 * length(r.amounts)
        @test alt["n_starts"][1] isa Int
        empty_sweep = SweepResult(sweep.case, [(T, P, "failed") for (T, P) in grid(sweep.case)])
        @test isempty(species_table(empty_sweep)) && csv_text(phase_summary(empty_sweep)) == "\n"
        @test csv_text(alternatives_table(empty_sweep)) ==
              "T_C,P_bar,result,gibbs_J_per_mol_C,n_starts,phase,species,mol_per_mol_C\n"
    end

    @testset "Python's files read as Python reads them" begin
        py = JSON.parsefile(joinpath(PY_RESULTS, "python_read.json"))
        tables = load_saved(PY_RUN, PY_RESULTS)
        @test tables.label == py["label"] && tables.case isa Case
        @test [[T, P, l] for ((T, P), l) in read_choices(PY_RUN)] == py["choices"]
        @test [[T, P] for (T, P) in flagged_points(tables)] == py["flagged_points"]
        @test same_as_python(tables.alternatives, py["alternatives"])
        for (k, name) in enumerate(("summary", "species", "phases"))
            @test same_as_python(getfield(tables, Symbol(name)), py["applied"][name])
            @test same_as_python(tables.base[k], py["base"][name])
        end
        # the chosen result B at 400 °C; Z (unknown) is ignored
        @test tables.choices == OrderedDict((300.0, 240.0) => "Z", (400.0, 240.0) => "B")
        @test tables.summary["result"] == ["A", "", "B", "A"]
        back = with_choices(tables, Dict())
        @test isequal(back.species, tables.base[2]) && back.summary["result"] == ["A", "", "A", "A"]
    end

    @testset "Julia's files read back to the same tables" begin
        dir = mktempdir()
        save(sweep, dir)
        for (f, table) in zip(RESULT_FILES, (summary(sweep), species_table(sweep), phase_summary(sweep),
                                            alternatives_table(sweep)))
            back = read_csv(joinpath(dir, f))
            @test back.columns == table.columns
            @test all(roundtrip_equal(back[c][i], table[c][i]) for c in table.columns for i in 1:nrows(table))
            @test csv_text(back) == read(joinpath(dir, f), String)   # read + write again: the same file
        end
        # the saved-run loader with choices written by Julia equals the one on Python's files
        write_choices(dir, Dict((400.0, 240.0) => "B", (300.0, 240.0) => "Z"))
        mine, theirs = load_saved(dir), load_saved(PY_RUN)
        for name in (:summary, :species, :phases, :alternatives)
            a, b = getfield(mine, name), getfield(theirs, name)
            @test a.columns == b.columns
            @test all(isequal(a[c][i], b[c][i]) || (a[c][i] isa Float64 && isapprox(a[c][i], b[c][i]; rtol=1e-14))
                      for c in a.columns for i in 1:nrows(a))
        end
        @test mine.choices == theirs.choices
    end
end
