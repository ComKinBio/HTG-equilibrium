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
                      resolve_user_dir
using Artifacts
using Dates
using Logging
using Sockets
using JSON
using Stipple, Stipple.ReactiveTools
using StippleUI
using Genie
using Electron

export State, worker_summary, status_text, output_folder, run_app, use_model!, ready

include("state.jl")      # State, worker_summary, output_folder, load!
include("page.jl")       # the page: header, status line, tabs; offline icon font
include("window.jl")     # can Electron open a window? (else the browser)
include("shell.jl")      # server, window, native dialogs, life cycle (run_app)
include("selftest.jl")   # HTG_SELFTEST=1: checks in the window
include("cli.jl")        # htg gui: main (@main)

end # module
