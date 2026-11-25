# Cross-check of the Case tab's files with the Python app: every case file of the
# app (cases/*.toml, the template, the old files in tests/data/legacy_cases) is
# opened in a CaseEditor, edited (a value changed and changed back) and saved
# with the tab's Save as; Python then loads each saved file and the original
# (htg.load_case) and compares the two cases as Python's dump_case writes them.
# A manual step, not part of Pkg.test() (it needs uv and the Python package):
#
#   julia --project=julia/HTGApp/test julia/HTGApp/test/crosscheck_python.jl [--out DIR]
#
# The saved files go to DIR (default: a new temporary folder); nothing is
# written to the repository. Not compared: the name (the saved files have other
# file names) and the [solver] keys fd_ftol and fd_eps, settings of the dropped
# solver slsqp-fd that the Julia app does not write (owner, PR #10). The
# comparison file is compared as the absolute path Python resolves.
using HTGApp, HTGEquilibrium
using HTGApp: edit!, write!, form_view

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))

const PYTHON = raw"""
import sys
from pathlib import Path
import htg
from htg.case import dump_case

def text(path):
    case = htg.load_case(path)
    if case.comparison:
        case.comparison = dict(case.comparison, file=str((Path(path).parent / case.comparison['file']).resolve()))
    # the name follows the file name when the file has none; the settings of the
    # dropped solver slsqp-fd are not written by the Julia app (owner, PR #10)
    return '\n'.join(line for line in dump_case(case).splitlines()
                     if not line.startswith(('name = ', 'fd_ftol = ', 'fd_eps = ')))

bad = 0
for line in sys.stdin.read().splitlines():
    original, saved = line.split('\t')
    a, b = text(original), text(saved)
    if a != b:
        bad += 1
        print(f'DIFFERENT: {original} -> {saved}')
        import difflib
        print('\n'.join(difflib.unified_diff(a.splitlines(), b.splitlines(), lineterm='')))
    else:
        print(f'same: {Path(original).relative_to(Path.cwd())}')
sys.exit(1 if bad else 0)
"""

function main(args)
    out = length(args) >= 2 && args[1] == "--out" ? mkpath(args[2]) : mktempdir()
    files = [joinpath(REPO, "cases", f) for f in readdir(joinpath(REPO, "cases")) if endswith(f, ".toml")]
    append!(files, [joinpath(REPO, "tests", "data", "legacy_cases", f)
                    for f in readdir(joinpath(REPO, "tests", "data", "legacy_cases")) if endswith(f, ".toml")])
    push!(files, HTGEquilibrium.CASE_TEMPLATE)
    pairs = String[]
    withenv("HTG_CONFIG_DIR" => mktempdir()) do
        state = State(; results_dir=out, cases_dir=joinpath(REPO, "cases"))
        use_model!(state, load_model(; user_data=nothing))   # an old species feed needs the molar masses
        for (i, file) in enumerate(files)
            ed = CaseEditor(state; async_bands=false)
            HTGApp.open!(ed, file) || error("could not open $file: $(ed.notifications)")
            if haskey(ed.fields, "solver.starts")
                starts = ed.texts["solver.starts"]
                edit!(ed, "solver.starts", "9")
                edit!(ed, "solver.starts", starts)
            end
            form_view(ed)["ok"] || error("$file: $(form_view(ed)["status"])")
            saved = joinpath(out, "$(i)_$(basename(file))")
            write!(ed, saved) || error("not saved: $(ed.notifications)")
            push!(pairs, "$file\t$saved")
        end
    end
    script = joinpath(out, "crosscheck.py")
    write(script, PYTHON)
    cmd = Cmd(`uv run python $script`; dir=REPO)
    ok = success(pipeline(cmd; stdin=IOBuffer(join(pairs, "\n")), stdout, stderr))
    println(ok ? "All $(length(pairs)) cases saved by the Case tab load in Python as the same case." :
                 "Some cases differ (above). Saved files: $out")
    return ok ? 0 : 1
end

exit(main(ARGS))
