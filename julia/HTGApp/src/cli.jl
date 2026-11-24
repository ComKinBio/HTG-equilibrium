# `htg gui` (Python: htg/cli.py cmd_gui). The launchers julia/htg and
# julia/htg.cmd run `julia --project=julia/HTGApp --threads=auto -m HTGApp ARGS...`
# when the command is `gui`. The arguments are parsed by the core's `parse_cli`
# (argparse's rules and messages, the same options as Python):
#
#     htg [--data DIR] gui [--browser] [--no-browser] [--workers N] [--port PORT] [--host HOST]
#                          [--results DIR] [--cases DIR] [--user-data DIR] [--no-user-species]
#
# Environment: HTG_SELFTEST=1 drives the window and checks it (selftest.jl);
# HTG_TIMING=1 prints time stamps of the start-up; HTG_T0 is the launch time as
# Unix seconds (date +%s.%N), so that Julia's own start-up counts too.

"""The `run_app` keywords for the parsed `gui` arguments (Python: cmd_gui)."""
function gui_options(values)
    user_data = resolve_user_dir(get(values, "user_data", nothing); enabled=!get(values, "no_user_species", false))
    return (results_dir=values["results"], cases_dir=values["cases"], workers=values["workers"],
            port=values["port"], show=!values["no_browser"], host=values["host"], browser=values["browser"],
            model_args=Dict{String,Any}("data_dir" => values["data"], "user_data" => user_data))
end

"""
    main(args=ARGS; out=stdout, err=stderr) -> exit code

`htg ... gui ...`: parse the arguments and run the app until its window is closed.
Other commands belong to the core (`julia -m HTGEquilibrium`).
"""
function main(args=ARGS; out::IO=stdout, err::IO=stderr)
    TIMING[] = get(ENV, "HTG_TIMING", "0") == "1" || get(ENV, "HTG_SELFTEST", "0") == "1"
    T0[] = something(tryparse(Float64, get(ENV, "HTG_T0", "")), time())
    stamp("julia running, HTGApp loaded")
    try
        command, values = HTGEquilibrium.parse_cli(collect(String, args), out, err)
        if command != "gui"
            println(err, "Error: HTGApp runs only the gui command; run the others with julia/htg $command")
            return 2
        end
        options = gui_options(values)
        on_usable = get(ENV, "HTG_SELFTEST", "0") == "1" ? selftest : nothing
        run_app(; options..., on_usable)
        return SELFTEST_FAILED[] ? 1 : 0
    catch exc
        exc isa HTGEquilibrium._CLIExit && return exc.code
        exc isa InterruptException && return 130
        showerror(err, exc, catch_backtrace())
        println(err)
        return 1
    end
end

@main
