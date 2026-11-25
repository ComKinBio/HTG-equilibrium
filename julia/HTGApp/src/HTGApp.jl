"""
HTGApp: the app of HTG equilibrium in Julia (Phase 2, milestone M4): a
Genie/Stipple page in an Electron window, on the science core HTGEquilibrium.
A separate package so that the core, its command line and its tests do not
load the GUI packages. Python counterpart: `src/htg/gui/`.

    julia/htg gui                                   # the launcher (Linux, macOS)
    julia --project=julia/HTGApp -t auto -m HTGApp gui [--browser] ...
"""
module HTGApp

using HTGEquilibrium
using HTGEquilibrium: Settings, load_settings, OrderedDict, load_model_from_args, worker_threads, all_molecular,
                      resolve_user_dir, feed_mode, water_molar_mass
using Artifacts
using Dates
using Printf
using Logging
using Sockets
using JSON
using Stipple, Stipple.ReactiveTools
using StippleUI
using Genie
using Electron

export State, worker_summary, status_text, output_folder, run_app, use_model!, ready
export CaseDocument, CaseEditor, open_document, new_document

include("case_form.jl")  # the Case tab's wording, field parsing, notices, derived numbers
include("document.jl")   # CaseDocument: the case being edited
include("state.jl")      # State, worker_summary, output_folder, load!
include("editor.jl")     # CaseEditor: the Case tab's form without the page
include("page.jl")       # the page: header, status line, tabs; offline icon font
include("case_tab.jl")   # the Case tab on the page (one CaseEditor per page)
include("window.jl")     # can Electron open a window? (else the browser)
include("shell.jl")      # server, window, native dialogs, life cycle (run_app)
include("selftest.jl")   # HTG_SELFTEST=1: checks in the window
include("cli.jl")        # htg gui: main (@main)

# Compile the Case tab's code into the package image. Otherwise the first page
# compiles it while the solver's warm-up compiles (about a minute), and the
# Case tab stays empty until both are done.
if ccall(:jl_generating_output, Cint, ()) == 1
    let dir = mktempdir()
        state = State(; results_dir=dir, cases_dir=dir, settings=HTGEquilibrium.Settings(joinpath(dir, "settings.toml")))
        ed = CaseEditor(state; async_bands=false, band_check=(T, P; margin_C=0.5) -> Any[])
        for (key, text) in (("feed.total", "10"), ("feed.solids", "15"), ("feed.element:C", "30.4"),
                            ("conditions.temperatures", "300:400 n=3"), ("aqueous.cutoff", "x"))
            edit!(ed, key, text)
        end
        share_unit!(ed, "ash", "kg")
        switch_conditions!(ed, "points")
        switch_feed!(ed, "species")
        JSON.json(form_view(ed))
        write!(ed, joinpath(dir, "case.toml"))
        case_tab_html()
        rm(dir; recursive=true, force=true)
    end
end

end # module
