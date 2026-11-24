# Server, Electron window, native dialogs, shut-down and memory report.
"A free TCP port on 127.0.0.1 (opened and closed again, so a tiny race is possible)."
function free_port()
    port, server = listenany(ip"127.0.0.1", 49152)
    close(server)
    Int(port)
end

"""
Run a JavaScript expression that gives a Promise in Electron's main process and wait
for its value. Electron.jl's `run(app, code)` only returns plain values, so the result
is parked in a global and polled.
"""
function electron_await(app, js; timeout = 3600)
    key = string(rand(UInt64))
    run(app, """(global.__m0 = global.__m0 || {}, global.__m0['$key'] = null,
                 Promise.resolve($js).then(r => global.__m0['$key'] = {ok: r},
                                           e => global.__m0['$key'] = {err: String(e)}), null)""")
    t = time()
    while time() - t < timeout
        r = run(app, "(() => { const r = global.__m0['$key']; if (r !== null) delete global.__m0['$key']; return r })()")
        if r !== nothing
            haskey(r, "err") && error("Electron: " * r["err"])
            return r["ok"]
        end
        sleep(0.05)
    end
    error("Electron did not answer within $timeout s")
end

function native_dialogs(app, win)
    parent = "electron.BrowserWindow.fromId($(win.id))"
    M0App.NATIVE_OPEN[] = opts -> begin
        r = electron_await(app, "electron.dialog.showOpenDialog($parent, $(JSON.json(opts)))")
        r["canceled"] || isempty(r["filePaths"]) ? nothing : String(r["filePaths"][1])
    end
    M0App.NATIVE_SAVE[] = opts -> begin
        r = electron_await(app, "electron.dialog.showSaveDialog($parent, $(JSON.json(opts)))")
        r["canceled"] || isempty(something(r["filePath"], "")) ? nothing : String(r["filePath"])
    end
end

# ---------------------------------------------------------------- memory (Linux)
function children_of(pid)
    kids = Int[]
    for d in readdir("/proc")
        all(isdigit, d) || continue
        stat = try read("/proc/$d/stat", String) catch; continue end
        fields = split(stat[findlast(')', stat)+2:end])
        parse(Int, fields[2]) == pid && push!(kids, parse(Int, d))
    end
    kids
end

function tree(pid)
    out = [pid]
    for k in children_of(pid)
        append!(out, tree(k))
    end
    out
end

function mem_kb(pid, key)
    for (file, k) in (("/proc/$pid/smaps_rollup", key),)
        isfile(file) || continue
        for line in eachline(file)
            startswith(line, k * ":") && return parse(Int, split(line)[2])
        end
    end
    0
end

function memory_report(app)
    jl = getpid()
    el = tree(getpid(app.proc))
    jl_rss, jl_pss = mem_kb(jl, "Rss"), mem_kb(jl, "Pss")
    el_rss, el_pss = sum(mem_kb(p, "Rss") for p in el), sum(mem_kb(p, "Pss") for p in el)
    mb(x) = round(Int, x / 1024)
    stamp("memory: julia RSS $(mb(jl_rss)) MB (PSS $(mb(jl_pss)) MB); " *
          "electron $(length(el)) processes RSS $(mb(el_rss)) MB (PSS $(mb(el_pss)) MB); " *
          "total PSS $(mb(jl_pss + el_pss)) MB")
end

# ---------------------------------------------------------------- screenshot

function screenshot(app, win, path)
    js = """electron.BrowserWindow.fromId($(win.id)).webContents.capturePage()
              .then(img => { require('fs').writeFileSync($(JSON.json(path)), img.toPNG()); return true })"""
    electron_await(app, js)
    stamp("screenshot written to $(path)")
end

# ---------------------------------------------------------------- main
function main(t0::Float64 = time())
    T0[] = t0
    ready = Channel{Float64}(1)
    M0App.ON_READY[] = _ -> (Base.isready(ready) || put!(ready, time()))

    port = free_port()
    Genie.config.server_host = "127.0.0.1"
    Genie.config.log_requests = false
    M0App.register()
    Genie.up(port, "127.0.0.1"; async = true, open_browser = false)
    stamp("server listening on 127.0.0.1:$port")

    # On GNOME, Electron's file dialogs go through the XDG desktop portal (GNOME's own
    # file chooser). M0_GTK_DIALOGS=1 forces Electron's in-process GTK dialog instead.
    args = get(ENV, "M0_GTK_DIALOGS", "0") == "1" ? ["--xdg-portal-required-version=99"] : String[]
    app = Electron.Application(; additional_electron_args = args)
    stamp("electron process started")
    opts = Dict("title" => "HTG equilibrium (M0)", "width" => 1400, "height" => 900,
                "webPreferences" => Dict("nodeIntegration" => false, "contextIsolation" => true))
    win = Window(app, URI("http://127.0.0.1:$(port)/"); options = opts)
    stamp("window shows the page (did-finish-load)")
    native_dialogs(app, win)

    timer = Timer(120)
    @async (wait(timer); Base.isready(ready) || put!(ready, NaN))
    t_ready = take!(ready)
    close(timer)
    isnan(t_ready) && error("the page did not connect within 120 s")
    println("[M0] $(lpad(round(t_ready - T0[]; digits = 2), 6)) s  USABLE: Vue mounted, websocket connected, server handler ran"); flush(stdout)

    get(ENV, "M0_SHOT", "") != "" && (sleep(1.5); screenshot(app, win, ENV["M0_SHOT"]))
    get(ENV, "M0_SELFTEST", "0") == "1" && selftest(app, win)
    idle = parse(Float64, get(ENV, "M0_IDLE", "0"))
    if idle > 0
        sleep(idle)
        memory_report(app)
    end
    # debugging aid: run every *.js file dropped into M0_JSDIR in the window, write *.out
    jsdir = get(ENV, "M0_JSDIR", "")
    if jsdir != ""
        @async while isopen(win)
            for f in filter(endswith(".js"), readdir(jsdir; join = true))
                code = read(f, String)
                rm(f)
                # name.main.js runs in Electron's main process, other files in the page
                target = endswith(f, ".main.js") ? app : win
                res = try JSON.json(run(target, code)) catch e sprint(showerror, e) end
                write(f[1:end-3] * ".out", res)
            end
            sleep(0.2)
        end
    end
    if get(ENV, "M0_AUTOCLOSE", "0") == "1"
        stamp("closing the window as a user would (BrowserWindow.close)")
        ElectronAPI.close(win)
    end

    # Electron.jl's main.js keeps Electron alive when the last window closes, so wait for
    # the window, then end Electron, the server and Julia ourselves.
    while isopen(win) && isopen(app)
        sleep(0.1)
    end
    stamp("window closed")
    isopen(app) && close(app)
    t = time()
    while process_running(app.proc) && time() - t < 10
        sleep(0.05)
    end
    process_running(app.proc) && kill(app.proc)
    Genie.Server.down!()
    stamp("electron exited ($(process_exited(app.proc))), server stopped")
    nothing
end

