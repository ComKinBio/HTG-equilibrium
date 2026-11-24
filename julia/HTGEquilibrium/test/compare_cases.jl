# Run the example cases (cases/*.toml) in both versions and compare the saved
# results: file names, columns, rows and values (migration plan: "Run all
# example cases in cases/ in both versions and compare summary.csv and
# species.csv"). A manual step, not part of Pkg.test() (it takes a few minutes):
#
#   JULIA_NUM_THREADS=4 julia --project=julia/HTGEquilibrium julia/HTGEquilibrium/test/compare_cases.jl \
#       [--out DIR] [--rerun-python] [--rerun-julia] [case.toml ...]
#
# DIR (default: a folder in the system's temporary directory) gets python/<case>
# and julia/<case>. Python runs with the command line, `uv run htg run CASE --out
# DIR/python/<case> --no-plots --no-user-species` (only when its folder does not
# exist yet, or with --rerun-python); Julia runs all cases on one WorkerPool
# (threads), which saves each run when it finishes, as the Python worker pool
# does. Generated outputs are not committed.
#
# Comparison rules (as the reference solve tests): strings exact; Gibbs energy
# per mol C rel 1e-9; amounts (mol per mol C) rtol 1e-3 / atol 1e-4; columns
# derived from an amount (per kg dry feed, g, carbon) with the same tolerance on
# the implied amount; carbon percentages rtol 1e-3 / atol 1e-2; phase Gibbs
# energies rel 1e-6 / atol 1 J/mol C. Solver diagnostics (residuals,
# iterations, descent_RT, carbon_bound_violation) are reported, not compared;
# residuals must be below 1e-7 in both. A point where Julia's G is lower than
# Python's is accepted when Julia's solution is feasible (balance residual) and
# stationary (descent_RT small); it is listed. case.toml must be equal apart
# from the [solver] keys fd_ftol/fd_eps (settings of the dropped solver
# slsqp-fd, still written by Python; owner, PR #10).
using HTGEquilibrium
using HTGEquilibrium: OrderedDict, py_format_g
using Printf

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))

function parse_args(args)
    out = joinpath(tempdir(), "htg_compare_cases")
    rerun_python = rerun_julia = false
    cases = String[]
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--out"
            out = args[i+1]; i += 1
        elseif a == "--rerun-python"
            rerun_python = true
        elseif a == "--rerun-julia"
            rerun_julia = true
        else
            push!(cases, abspath(a))
        end
        i += 1
    end
    isempty(cases) && (cases = sort([joinpath(REPO, "cases", f) for f in readdir(joinpath(REPO, "cases"))
                                     if endswith(f, ".toml")]))
    return abspath(out), rerun_python, rerun_julia, cases
end

function run_python(case_file, dir)
    log = dir * ".log"
    cmd = `uv run htg run $case_file --out $dir --no-plots --no-user-species`
    t = @elapsed ok = success(pipeline(Cmd(cmd; dir=REPO); stdout=log, stderr=log))
    println("  Python: $(basename(case_file)) in $(round(t; digits=1)) s (log: $log)")
    ok || println("  Python exit code != 0 (failed points; see the log)")
    return t
end

function run_julia(cases, dirs)
    model = load_model()
    t = @elapsed begin
        WorkerPool(min(Threads.nthreads(), 4); model, warm_up=false) do pool
            runs = [submit(pool, load_case(c), d) for (c, d) in zip(cases, dirs)]
            foreach(wait, runs)
            for r in runs
                last(r.log) == "results saved to $(r.output_dir)" || println("  Julia: ", r.name, ": ", last(r.log))
            end
        end
    end
    println("  Julia: $(length(cases)) cases on $(Threads.nthreads()) threads in $(round(t; digits=1)) s " *
            "(including compilation)")
    return t
end

# ---------------------------------------------------------------- comparison
isnum(x) = x isa Real && !(x isa Bool)
fl(x) = isnum(x) ? Float64(x) : NaN
nullish(x) = x === missing || (x isa AbstractFloat && isnan(x)) || x == ""
close_to(a, b; rtol, atol) = (isnan(a) && isnan(b)) || abs(a - b) <= atol + rtol * abs(b)

mutable struct Report
    points::Int
    matched::Int
    lower_G::Vector{String}
    problems::Vector{String}
    worst::OrderedDict{String,Float64}   # column => worst |a-b| / (atol + rtol|b|)
    notes::Vector{String}
end
Report() = Report(0, 0, String[], String[], OrderedDict{String,Float64}(), String[])

function note_worst!(rep, key, a, b; rtol, atol)
    (isnan(a) && isnan(b)) && return true
    ratio = abs(a - b) / (atol + rtol * abs(b))
    rep.worst[key] = max(get(rep.worst, key, 0.0), isnan(ratio) ? Inf : ratio)
    return ratio <= 1
end

const AMOUNT = (rtol=1e-3, atol=1e-4)
const PERCENT = (rtol=1e-3, atol=1e-2)
const GIBBS = (rtol=1e-9, atol=0.0)
const PHASE_GIBBS = (rtol=1e-6, atol=1.0)
const DIAGNOSTIC = ("balance_residual", "charge_residual", "iterations", "descent_RT", "carbon_bound_violation")
const DERIVED = ("mol_per_kg_dry_feed", "g_per_kg_dry_feed", "carbon_mol_per_kg_dry_feed")

pointkey(t, i) = (fl(t["T_C"][i]), haskey(t, "P_bar") ? fl(t["P_bar"][i]) : NaN)
label(p) = "$(py_format_g(p[1], 6)) °C, $(py_format_g(p[2], 6)) bar"

"""Compare one table row by row; returns the set of points with a difference."""
function compare_rows!(rep, name, a, b, bad_points; gibbs_lower=Set{Tuple{Float64,Float64}}())
    for i in 1:nrows(b)
        p = pointkey(b, i)
        p == pointkey(a, i) || (push!(rep.problems, "$name row $i: point $(label(pointkey(a, i))) vs $(label(p))");
                                push!(bad_points, p); continue)
        p in gibbs_lower && continue   # another minimum: compared by the lower-G check
        for c in b.columns
            x, y = a[c][i], b[c][i]
            c in DIAGNOSTIC && continue
            ok = true
            if !(isnum(x) && isnum(y))
                ok = isequal(x, y) || (nullish(x) && nullish(y)) || (isnum(x) && isnum(y) && fl(x) == fl(y))
                ok || push!(rep.problems, "$name $(label(p)) $c: Julia $(repr(x)) vs Python $(repr(y))")
            elseif c in ("T_C", "P_bar")
                ok = x == y
            elseif c == "gibbs_J_per_mol_C" && name != "phases"
                ok = note_worst!(rep, "$name.$c", fl(x), fl(y); GIBBS...)
            elseif c == "gibbs_J_per_mol_C" || c == "phase_gibbs_J_per_mol_C"
                ok = note_worst!(rep, "$name.$c", fl(x), fl(y); PHASE_GIBBS...)
            elseif c in ("carbon_percent", "gas_carbon_percent")
                ok = note_worst!(rep, "$name.$c", fl(x), fl(y); PERCENT...)
            elseif c in DERIVED && haskey(b, "mol_per_mol_C")
                mol = fl(b["mol_per_mol_C"][i])
                if mol != 0 && fl(y) != 0
                    k = fl(y) / mol   # implied-amount scale of this row
                    ok = note_worst!(rep, "$name.$c", fl(x) / k, fl(y) / k; AMOUNT...)
                else
                    ok = note_worst!(rep, "$name.$c", fl(x), fl(y); rtol=AMOUNT.rtol, atol=1e-2)
                end
            else
                ok = note_worst!(rep, "$name.$c", fl(x), fl(y); AMOUNT...)
            end
            if !ok
                push!(bad_points, p)
                push!(rep.problems, "$name $(label(p)) $c: Julia $x vs Python $y")
            end
        end
    end
end

function compare_case(case_name, jdir, pdir)
    rep = Report()
    jfiles, pfiles = sort(readdir(jdir)), sort(readdir(pdir))
    jfiles == pfiles || push!(rep.problems, "files: Julia $jfiles vs Python $pfiles")
    # case.toml: equal except the dropped solver's keys, which Python still writes
    ptext = join([l for l in split(read(joinpath(pdir, "case.toml"), String), '\n')
                  if !startswith(l, "fd_ftol") && !startswith(l, "fd_eps")], '\n')
    ptext == read(joinpath(jdir, "case.toml"), String) || push!(rep.problems, "case.toml differs")
    occursin("fd_", read(joinpath(pdir, "case.toml"), String)) &&
        push!(rep.notes, "case.toml: Python writes [solver] fd_ftol/fd_eps (dropped solver slsqp-fd), Julia does not")
    if isfile(joinpath(pdir, "comparison.csv"))
        read(joinpath(pdir, "comparison.csv")) == read(joinpath(jdir, "comparison.csv")) ||
            push!(rep.problems, "comparison.csv differs")
    end
    tables = OrderedDict(f => (read_csv(joinpath(jdir, f)), read_csv(joinpath(pdir, f)))
                         for f in ("summary.csv", "species.csv", "phases.csv", "alternatives.csv"))
    for (f, (a, b)) in tables
        a.columns == b.columns || push!(rep.problems, "$f columns: Julia $(a.columns) vs Python $(b.columns)")
        nrows(a) == nrows(b) || push!(rep.problems, "$f rows: Julia $(nrows(a)) vs Python $(nrows(b))")
    end
    !isempty(rep.problems) && return rep
    sa, sb = tables["summary.csv"]
    rep.points = nrows(sb)
    bad = Set{Tuple{Float64,Float64}}()
    lower = Set{Tuple{Float64,Float64}}()
    for i in 1:nrows(sb)
        p = pointkey(sb, i)
        sb["status"][i] == "solved" == sa["status"][i] || continue
        for c in ("balance_residual",)
            (abs(fl(sa[c][i])) <= 1e-7 && abs(fl(sb[c][i])) <= 1e-7) ||
                push!(rep.problems, "summary $(label(p)) $c: Julia $(sa[c][i]), Python $(sb[c][i])")
        end
        Gj, Gp = fl(sa["gibbs_J_per_mol_C"][i]), fl(sb["gibbs_J_per_mol_C"][i])
        if Gj < Gp - 1e-9 * abs(Gp)
            stationary = fl(sa["descent_RT"][i]) <= 1e-3 && abs(fl(sa["balance_residual"][i])) <= 1e-7
            push!(rep.lower_G, @sprintf("%s: Julia G lower by %.3g (rel %.2g); Julia descent_RT %.2g, residual %.2g%s",
                                        label(p), Gp - Gj, (Gp - Gj) / abs(Gp), fl(sa["descent_RT"][i]),
                                        fl(sa["balance_residual"][i]), stationary ? " (verified)" : " (NOT stationary)"))
            push!(lower, p)
            stationary || push!(bad, p)
        end
        sa["iterations"][i] == sb["iterations"][i] ||
            push!(rep.notes, "$(label(p)): iterations Julia $(sa["iterations"][i]), Python $(sb["iterations"][i])")
    end
    compare_rows!(rep, "summary", sa, sb, bad; gibbs_lower=lower)
    for (f, name) in (("species.csv", "species"), ("phases.csv", "phases"), ("alternatives.csv", "alternatives"))
        compare_rows!(rep, name, tables[f]..., bad; gibbs_lower=lower)
    end
    rep.matched = rep.points - length(union(bad, lower))
    return rep
end

function main(args)
    out, rerun_python, rerun_julia, cases = parse_args(args)
    names = [load_case(c).name for c in cases]
    pdirs = [joinpath(out, "python", n) for n in names]
    jdirs = [joinpath(out, "julia", n) for n in names]
    println("Output: $out")
    for (c, d) in zip(cases, pdirs)
        (rerun_python || !isdir(d)) && run_python(c, d)
    end
    todo = [i for i in eachindex(cases) if rerun_julia || !isdir(jdirs[i])]
    isempty(todo) || run_julia(cases[todo], jdirs[todo])
    all_ok = true
    for (n, j, p) in zip(names, jdirs, pdirs)
        rep = compare_case(n, j, p)
        ok = isempty(rep.problems)
        all_ok &= ok
        println("\n$n: $(rep.points) points, $(rep.matched) match, $(length(rep.lower_G)) with a lower G in Julia, " *
                "$(length(rep.problems)) problems")
        for line in rep.lower_G
            println("  lower G: ", line)
        end
        for line in first(rep.problems, 30)
            println("  PROBLEM: ", line)
        end
        length(rep.problems) > 30 && println("  ... $(length(rep.problems) - 30) more")
        for line in rep.notes
            println("  note: ", line)
        end
        worst = sort(collect(rep.worst); by=last, rev=true)
        println("  worst deviation / tolerance: ",
                join(["$k $(@sprintf("%.2g", v))" for (k, v) in first(worst, 6)], ", "))
    end
    println(all_ok ? "\nAll cases agree." : "\nDifferences found (see above).")
    return all_ok
end

main(ARGS) || exit(1)
