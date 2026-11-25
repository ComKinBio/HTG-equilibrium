# Tests of HTGApp that need no display: the status line text, State, settings,
# the gui arguments, the window check, and the page as the server sends it.
#
#   julia --project=julia/HTGApp -t 4 -e 'using Pkg; Pkg.test()'
#
# The window itself is checked by the self-test (test/selftest.jl, needs a desktop).
using HTGApp, HTGEquilibrium, Test
using HTGApp: gui_options, explain_electron_failure, window_or_browser, resolve_path, placeholder_text, TABS,
              load!, STATE, register_page, FONT_ROUTE, icon_font_path, local_url
using HTGEquilibrium: OrderedDict, JSON, parse_cli, load_settings
import HTGApp.Genie, HTGApp.Dates

const PY_STATE = JSON.parsefile(joinpath(@__DIR__, "data", "state_python.json"))

"""Run HTGApp.main(args) and return (exit code, stdout, stderr)."""
function run_main(args)
    out, err = IOBuffer(), IOBuffer()
    code = HTGApp.main(args; out, err)
    return code, String(take!(out)), String(take!(err))
end

"""A pool with a stand-in solver (no equation of state), as the core's queue tests."""
stand_in(model, case, chain, on_point; stop=() -> false) = (sleep(0.05); for (T, P) in chain; on_point(T, P, "stand-in", ""); end)

@testset "HTGApp" begin
    @testset "worker_summary (Python's strings)" begin
        for c in PY_STATE["worker_summary"]
            status = OrderedDict{String,Any}("workers" => [OrderedDict{String,Any}("id" => i, "state" => s)
                                                           for (i, s) in enumerate(c["states"])],
                                             "queued_pieces" => c["queued_pieces"], "target" => length(c["states"]))
            @test worker_summary(status) == c["text"]
        end
    end

    @testset "output folder" begin
        @test output_folder("/r", "my case", Dates.DateTime(2026, 11, 4, 9, 5, 7)) == joinpath("/r", "my case-20261104-090507")
        @test startswith(basename(output_folder("results", "x")), "x-20")
    end

    @testset "paths resolved as Python's Path.resolve" begin
        dir = realpath(mktempdir())
        @test resolve_path(joinpath(dir, "a", "..", "b")) == joinpath(dir, "b")
        cd(dir) do
            @test resolve_path("results") == joinpath(dir, "results")
        end
        symlink(dir, joinpath(dir, "link"))
        @test resolve_path(joinpath(dir, "link", "new")) == joinpath(dir, "new")
    end

    config = mktempdir()
    withenv("HTG_CONFIG_DIR" => config) do
        @testset "State" begin
            work = mktempdir()
            state = State(; results_dir=joinpath(work, "res"), cases_dir=joinpath(work, "cases"))
            @test isabspath(state.results_dir) && basename(state.results_dir) == "res" && basename(state.cases_dir) == "cases"
            @test state.settings.path == joinpath(config, "settings.toml")   # Settings.load() by default
            @test !ready(state) && status_text(state) == "Loading the species data"
            @test_throws ArgumentError submit(state, load_case(HTGEquilibrium.CASE_TEMPLATE))
            # settings: read once, written by the tabs
            state.settings["temperature_unit"] = "K"
            @test load_settings()["temperature_unit"] == "K"
            @test State(; settings=load_settings(mktempdir())).settings["temperature_unit"] == "C"

            model = load_model(; user_data=nothing)
            use_model!(state, model)
            @test state.model === model
            @test Dict(k => v for (k, v) in state.phase_options) == PY_STATE["phase_options"]
            @test collect(keys(state.phase_options)) == ["gas", "aqueous", "organic_liquid", "solid"]
            @test state.feed_species == PY_STATE["feed_species"]

            # with a pool: the status line is worker_summary; submit puts the run below results_dir
            pool = WorkerPool(2; model, warm_up=false, chain_solver=stand_in)
            state.pool = pool
            HTGEquilibrium.start(pool)
            @test ready(state)
            timedwait(() -> status_text(state) == "Workers: 2 idle", 120.0)   # the first call compiles
            @test status_text(state) == "Workers: 2 idle"
            case = load_case(HTGEquilibrium.CASE_TEMPLATE)
            run = submit(state, case)
            @test startswith(run.output_dir, joinpath(state.results_dir, case.name * "-"))
            @test wait(run, 30)
            shutdown(pool)
        end

        @testset "loading in the background" begin
            state = State()
            out = IOBuffer()
            load!(state; workers=1, warm_up=false, model_args=Dict("user_data" => nothing), out)
            @test ready(state) && state.pool.target == 1 && isempty(state.problem)
            @test occursin("Starting 1 worker thread(s)", String(take!(out)))
            @test status_text(state) == "Workers: 1 idle" || timedwait(() -> status_text(state) == "Workers: 1 idle", 60.0) === :ok
            shutdown(state.pool)
            # a broken data folder: the status line says so
            bad = State()
            redirect_stderr(devnull) do
                load!(bad; model_args=Dict("data_dir" => joinpath(config, "missing"), "user_data" => nothing), out)
            end
            @test !ready(bad) && !isempty(bad.problem)
            @test startswith(status_text(bad), "Could not load the species data: ")
        end
    end

    @testset "gui arguments (Python's htg gui options)" begin
        _, v = parse_cli(["gui"])
        o = gui_options(v)
        @test o.results_dir == "results" && o.cases_dir == "cases" && o.workers === nothing && o.port === nothing
        @test o.show && o.host == "127.0.0.1" && !o.browser && o.model_args["data_dir"] === nothing
        @test o.model_args["user_data"] == HTGEquilibrium.default_user_dir()
        _, v = parse_cli(["--data", "d", "gui", "--browser", "--no-b", "--workers", "3", "--port=9000", "--host", "0.0.0.0",
                          "--results", "r", "--cases", "c", "--no-user-species"])
        o = gui_options(v)
        @test o.browser && !o.show && o.workers == 3 && o.port == 9000 && o.host == "0.0.0.0"
        @test o.results_dir == "r" && o.cases_dir == "c" && o.model_args == Dict("data_dir" => "d", "user_data" => nothing)
        @test gui_options(parse_cli(["gui", "--user-data", "~/u"])[2]).model_args["user_data"] == expanduser("~/u")
        # usage errors and help as the core's command line (argparse)
        code, out, err = run_main(["gui", "--port", "x"])
        @test code == 2 && endswith(err, "htg gui: error: argument --port: invalid int value: 'x'\n")
        code, out, err = run_main(["gui", "--help"])
        @test code == 0 && startswith(out, "usage: htg gui [-h] [--browser]")
        code, out, err = run_main(["species"])
        @test code == 2 && startswith(err, "Error: HTGApp runs only the gui command")
        @test local_url("0.0.0.0", 8080) == "http://127.0.0.1:8080/" && local_url("127.0.0.1", 5) == "http://127.0.0.1:5/"
    end

    @testset "window or browser (Python: window.py)" begin
        @test explain_electron_failure("[1:ERROR:ozone_platform_x11.cc(246)] Missing X server or \$DISPLAY\n" *
                                       "The platform failed to initialize.  Exiting.") ==
              "no display was found. Run htg gui from a desktop session."
        @test startswith(explain_electron_failure("electron: error while loading shared libraries: libnss3.so: " *
                                                  "cannot open shared object file: No such file or directory"),
                         "the system library libnss3.so is missing. Install it with your package manager")
        @test explain_electron_failure("x\nsomething else broke\n\n") == "Electron could not start (something else broke)."
        @test explain_electron_failure("") == "Electron could not start (Electron stopped without a message)."
        out = IOBuffer()
        @test window_or_browser(true; out, problem=() -> error("not asked")) == false && isempty(take!(out))
        @test window_or_browser(false; out, problem=() -> nothing) == true && isempty(take!(out))
        @test window_or_browser(false; out, problem=() -> "no display was found. Run htg gui from a desktop session.") == false
        @test String(take!(out)) == "Cannot open a separate window: no display was found. Run htg gui from a desktop " *
                                    "session. Opening the app in the browser instead (htg gui --browser does this directly).\n"
        withenv("DISPLAY" => nothing, "WAYLAND_DISPLAY" => nothing) do
            Sys.islinux() && @test HTGApp.window_problem() == "no display was found. Run htg gui from a desktop session."
        end
    end

    @testset "page (served offline, no .theme file)" begin
        @test [t[2] for t in TABS] == ["Case", "Queue", "Results", "Compare", "Species"]
        @test [t[3] for t in TABS] == ["edit_note", "queue", "insights", "compare_arrows", "science"]
        @test placeholder_text("Queue") == "The Queue tab is not in the Julia version yet. For now, use the Python app (uv run htg gui)."
        @test filesize(icon_font_path()) > 100_000
        dir = mktempdir()
        cd(dir) do
            STATE[] = State(; settings=load_settings(mktempdir()))
            Genie.config.app_env = "prod"
            Genie.config.log_requests = false
            HTGApp.private_sessions!()
            register_page()
            port = HTGApp.free_port()
            Genie.up(port, "127.0.0.1"; async=true, open_browser=false)
            try
                page = String(Genie.Router.HTTP.get("http://127.0.0.1:$port/").body)
                @test occursin("<title>HTG equilibrium</title>", page)
                @test !occursin("googleapis", page) && !occursin("fonts.g", page) && !occursin("/css/themes/", page)
                @test !occursin(r"(href|src)=['\"]https?://", page)    # nothing loaded from the internet
                @test occursin("quasar.prod", page) && occursin(FONT_ROUTE, page) && !occursin("stipplecore.css", page)
                for (_, label, icon) in TABS
                    @test occursin("label=\"$label\"", page) && occursin("icon=\"$icon\"", page)
                end
                font = Genie.Router.HTTP.get("http://127.0.0.1:$port$FONT_ROUTE")
                @test font.status == 200 && font.body == read(icon_font_path())
                # the session cookie: a name of this run, no expiry date; a cookie of an earlier run is ignored
                sessions = HTGApp.Stipple.ModelStorage.Sessions
                name = sessions.GenieSession.session_key_name()
                @test startswith(name, "__htg_") && !haskey(sessions.GenieSession.session_options(), "MaxAge")
                @test startswith(basename(sessions.GenieSessionFileSession.sessions_path()), "htg-sessions-")
                response = Genie.Router.HTTP.get("http://127.0.0.1:$port/")
                cookie = join([v for (k, v) in response.headers if lowercase(k) == "set-cookie"], "; ")
                @test occursin(name * "=", cookie) && !occursin(r"max-age"i, cookie)
                logs, _ = Test.collect_test_logs() do
                    Genie.Router.HTTP.get("http://127.0.0.1:$port/"; headers=["Cookie" => "__geniesid=$("ab"^40); __htg_0123456789abcdef=$("cd"^40)"])
                end
                @test isempty(logs)
            finally
                Genie.Server.down!()
            end
            @test isempty(readdir(dir))   # no .theme (Stipple's theme dotfile) or anything else
        end
    end
end
