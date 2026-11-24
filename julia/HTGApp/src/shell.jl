# Server, window and life cycle (Python: gui/app.py run_app, gui/window.py).
# Taken from the M0 spike (julia/spikes/m0_window/M0App/src/shell.jl).

"""Launch time (Unix seconds) that the time stamps are relative to (`HTG_T0`, else the start of `run_app`)."""
const T0 = Ref(time())
"""Print time stamps of the start-up (`HTG_TIMING=1`; the self-test turns it on)."""
const TIMING = Ref(false)
stamp(msg) = TIMING[] && (println("[htg] $(lpad(round(time() - T0[]; digits = 2), 6)) s  ", msg); flush(stdout))

"""A free TCP port on 127.0.0.1 (opened and closed again, so a tiny race is possible)."""
function free_port()
    port, server = listenany(ip"127.0.0.1", 49152)
    close(server)
    return Int(port)
end

"""The address to open: 0.0.0.0 and :: listen everywhere, this computer is 127.0.0.1."""
local_url(host, port) = "http://$(host in ("0.0.0.0", "::", "") ? "127.0.0.1" : host):$port/"

"""Open `url` in the default web browser (no error if that fails)."""
function open_in_browser(url)
    cmd = Sys.iswindows() ? `cmd /c start "" $url` : Sys.isapple() ? `open $url` : `xdg-open $url`
    try
        run(pipeline(ignorestatus(cmd); stdout=devnull, stderr=devnull); wait=false)
    catch
    end
    return nothing
end

# ---------------------------------------------------------------- native dialogs
"""
Run a JavaScript expression that gives a Promise in Electron's main process and wait
for its value. Electron.jl's `run(app, code)` only returns plain values, so the result
is parked in a global and polled (M0 spike).
"""
function electron_await(app, js; timeout = 3600)
    key = string(rand(UInt64))
    run(app, """(global.__htg = global.__htg || {}, global.__htg['$key'] = null,
                 Promise.resolve($js).then(r => global.__htg['$key'] = {ok: r},
                                           e => global.__htg['$key'] = {err: String(e)}), null)""")
    t = time()
    while time() - t < timeout
        r = run(app, "(() => { const r = global.__htg['$key']; if (r !== null) delete global.__htg['$key']; return r })()")
        if r !== nothing
            haskey(r, "err") && error("Electron: " * r["err"])
            return r["ok"]
        end
        sleep(0.05)
    end
    error("Electron did not answer within $timeout s")
end

"""
Native file dialogs of the window, for the tabs (owner decision, PR #3): each takes
Electron's dialog options as a Dict and returns a path or `nothing` (cancelled, or
no window: browser mode).
"""
const NATIVE_OPEN = Ref{Any}(nothing)
const NATIVE_SAVE = Ref{Any}(nothing)
open_dialog(options) = NATIVE_OPEN[] === nothing ? nothing : NATIVE_OPEN[](options)
save_dialog(options) = NATIVE_SAVE[] === nothing ? nothing : NATIVE_SAVE[](options)

function native_dialogs!(app, win)
    parent = "electron.BrowserWindow.fromId($(win.id))"
    NATIVE_OPEN[] = opts -> begin
        r = electron_await(app, "electron.dialog.showOpenDialog($parent, $(JSON.json(opts)))")
        r["canceled"] || isempty(r["filePaths"]) ? nothing : String(r["filePaths"][1])
    end
    NATIVE_SAVE[] = opts -> begin
        r = electron_await(app, "electron.dialog.showSaveDialog($parent, $(JSON.json(opts)))")
        r["canceled"] || isempty(something(r["filePath"], "")) ? nothing : String(r["filePath"])
    end
    return nothing
end

"""
Let the window load only this server's pages: every other request (fonts, scripts
from the internet, links) is cancelled and remembered in `global.__htgBlocked`
(the self-test checks that it stays empty).
"""
function local_only!(app, port)
    run(app, """(() => {
      global.__htgBlocked = [];
      electron.session.defaultSession.webRequest.onBeforeRequest((details, callback) => {
        let ok = false;
        try {
          const u = new URL(details.url);
          ok = ['data:', 'devtools:', 'blob:'].includes(u.protocol) ||
               (['http:', 'ws:'].includes(u.protocol) && u.hostname === '127.0.0.1' && u.port === '$port');
        } catch (e) {}
        if (!ok) global.__htgBlocked.push(details.url);
        callback({cancel: !ok});
      });
      return null; })()""")
end

blocked_requests(app) = run(app, "global.__htgBlocked || []")

# Shown at once while the server compiles the real page (without a sysimage that takes ~40 s).
const LOADING_HTML = """
<!DOCTYPE html><html><head><meta charset="utf-8"><title>HTG equilibrium</title><style>
body { margin: 0; font-family: Roboto, -apple-system, "Helvetica Neue", Helvetica, Arial, sans-serif; }
header { background: #5898d4; color: white; padding: 12px 16px 14px; font-size: 1.25rem; font-weight: 500;
  box-shadow: 0 1px 5px rgba(0,0,0,.2); }
p { margin: 24px 16px; font-size: 1rem; color: #555; }
</style></head><body><header>HTG equilibrium</header><p>Starting the app.</p></body></html>"""

# ---------------------------------------------------------------- the window
"""
    run_window(url, port; on_usable=nothing) -> Bool

Open the Electron window (1400 x 900, "HTG equilibrium") on `url` and wait until
the user closes it; false (after one printed line) when Electron cannot start.
`on_usable(app, win)` runs once the page is connected (the self-test).
"""
function run_window(url, port; on_usable=nothing)
    ready = Channel{Float64}(1)
    ON_READY[] = _ -> (Base.isready(ready) || put!(ready, time()))
    app = try
        Electron.Application()
    catch exc
        exc isa InterruptException && rethrow()
        first_line = first(split(sprint(showerror, exc), '\n'))
        println("Cannot open a separate window: Electron could not start ($first_line). " *
                "Opening the app in the browser instead (htg gui --browser does this directly).")
        return false
    end
    stamp("electron process started")
    local_only!(app, port)
    opts = Dict("title" => "HTG equilibrium", "width" => 1400, "height" => 900,
                "webPreferences" => Dict("nodeIntegration" => false, "contextIsolation" => true))
    win = Window(app, LOADING_HTML; options = opts)
    stamp("window open (start-up page)")
    native_dialogs!(app, win)
    @async try
        load(win, URI(url))
        stamp("window shows the app page (did-finish-load)")
    catch exc
        isopen(win) && @warn "Could not load the app page" exception = exc
    end
    @async begin
        t_ready = take!(ready)
        println("[htg] $(lpad(round(t_ready - T0[]; digits = 2), 6)) s  USABLE: page loaded, websocket connected, server handler ran")
        flush(stdout)
        on_usable === nothing || on_usable(app, win)
    end
    # Electron.jl's main.js keeps Electron alive when the last window closes, so wait for
    # the window, then end Electron ourselves.
    try
        while isopen(win) && isopen(app)
            sleep(0.1)
        end
    catch exc
        exc isa InterruptException || rethrow()   # Ctrl+C in the terminal: close the window
    end
    stamp("window closed")
    NATIVE_OPEN[] = NATIVE_SAVE[] = nothing
    try
        isopen(app) && close(app)
    catch
    end
    t = time()
    while process_running(app.proc) && time() - t < 10
        sleep(0.05)
    end
    process_running(app.proc) && kill(app.proc)
    stamp("electron exited")
    return true
end

"""Serve until Ctrl+C (browser mode)."""
function serve_until_interrupted()
    Base.exit_on_sigint(false)
    try
        while true
            sleep(1)
        end
    catch exc
        exc isa InterruptException || rethrow()
    end
    return nothing
end

const SESSIONS_DIR = Ref("")

"""
Genie's session cookie and session files, for this run only. Genie encrypts the
cookie with the secret token, which is new in every run (`register_page`), and
keeps the cookie 30 days; a cookie of an earlier run then decrypts to garbage
and every request logs "Failed to store session data". So the cookie gets a new
name per run and lives only as long as the browser (or window). The session files
go to a temporary folder that is removed when Julia exits (Genie's default is a
fixed folder in the system temp folder that is never emptied).
"""
function private_sessions!()
    sessions = Stipple.ModelStorage.Sessions
    SESSIONS_DIR[] = mktempdir(; prefix = "htg-sessions-")
    sessions.GenieSessionFileSession.sessions_path(SESSIONS_DIR[])
    sessions.GenieSession.session_key_name("__htg_" * bytes2hex(rand(UInt8, 8)))
    options = copy(sessions.GenieSession.session_options())
    delete!(options, "MaxAge")
    sessions.GenieSession.session_options(options)
    return nothing
end

# ---------------------------------------------------------------- exit
"""
End the process now with `code`, without Julia's exit. Called by `main` after
`run_app` (the window, the server and the pool are stopped and "Stopped. ..."
is printed). Julia's exit waits until a compilation step in progress on a
worker thread has finished; without a sysimage that took 20-30 s when the
window was closed during the warm-up (and any first call of uncompiled code
waits in the same way). The session folder, the only temporary folder left,
is removed first.
"""
function exit_now(code)
    stamp("exiting")
    isempty(SESSIONS_DIR[]) || rm(SESSIONS_DIR[]; recursive = true, force = true)
    flush(stdout)
    flush(stderr)
    ccall(:_exit, Cvoid, (Cint,), code)
end

# ---------------------------------------------------------------- run_app
"""
    run_app(; results_dir="results", cases_dir="cases", workers=nothing, port=nothing, show=true,
            model_args=Dict(), host="127.0.0.1", browser=false, settings=nothing, warm_up=true,
            on_usable=nothing)

Open the app and serve it until the window is closed (or Ctrl+C in browser mode);
then stop the worker pool and the server (Python: `run_app`).

- `browser`: false (default): the app's own window (Electron); if that cannot
  start, print why and use the browser. true: a browser tab.
- `show`: browser mode only: open a browser tab (false: only serve).
- `port`: default 8080 in the browser, a free port for the window.
- `host`: "127.0.0.1" (default): only this computer can open the app. The app
  reads and writes files, so share it on a network only on purpose.
- `workers`: worker threads (default `default_workers()`).
- `model_args`: `load_model` arguments ("data_dir", "user_data", ...).

The window opens first; the model is loaded and the pool started in a
background task, and the status line says when it is ready.
"""
function run_app(; results_dir="results", cases_dir="cases", workers=nothing, port=nothing, show=true,
                 model_args=Dict{String,Any}(), host="127.0.0.1", browser=false, settings=nothing, warm_up=true,
                 on_usable=nothing)
    # Warnings and errors only: Genie's and Stipple's start-up messages are not for the user.
    global_logger(ConsoleLogger(stderr, Logging.Warn))
    Base.exit_on_sigint(false)   # Ctrl+C in the terminal stops the app the same way as closing the window
    # The window check starts a small Electron program; it runs while the page is registered.
    problem = browser || !show ? nothing : @async window_problem()
    state = State(; results_dir, cases_dir, settings)
    STATE[] = state
    # Compiled now: called during the warm-up, it would wait for the compiling worker.
    precompile(exit_now, (Int,))
    Threads.@spawn :default load!(state; workers, model_args, warm_up)
    stopped = Ref(false)
    function stop()
        stopped[] && return
        stopped[] = true
        state.pool === nothing || shutdown(state.pool; wait=false)
        try
            Genie.Server.down!()
        catch
        end
        println("Stopped. Finished runs are saved in ", state.results_dir)
        flush(stdout)
    end
    try
        ENV["GENIE_ENV"] = Genie.config.app_env = "prod"   # minified Vue/Quasar, no file watching
        Genie.config.log_requests = false
        Genie.config.server_host = host
        private_sessions!()
        register_page()
        window = window_or_browser(browser || !show; problem = () -> fetch(problem))
        stamp("page registered, window check done")
        port = something(port, window ? free_port() : 8080)
        # Genie prints empty lines around its (hidden) start-up messages.
        redirect_stdout(devnull) do
            Genie.up(port, host; async = true, open_browser = false)
        end
        url = local_url(host, port)
        stamp("server listening on $host:$port")
        if window && run_window(url, port; on_usable)
            return nothing
        end
        println("HTG equilibrium is running at $url (Ctrl+C stops it).")
        flush(stdout)
        show && open_in_browser(url)
        serve_until_interrupted()
    finally
        stop()
    end
    return nothing
end
