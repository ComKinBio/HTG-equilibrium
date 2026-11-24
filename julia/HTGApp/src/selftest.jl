# Self-test in the Electron window (HTG_SELFTEST=1; needs a desktop session).
# Drives the page with JavaScript, as the M0 spike did, prints one PASS/FAIL line
# per check, then closes the window as a user would. test/selftest.jl runs it
# twice (a setting written in the first session must be read in the second) and
# checks that no process is left.
#
#   HTG_SELFTEST_SET=unit      set the pressure unit preference to `unit` in this session
#   HTG_SELFTEST_EXPECT=unit   the pressure unit this session must have read at start-up
#   HTG_SELFTEST_EARLY=1       close the window as soon as the page is checked (during the warm-up)
#   HTG_SELFTEST_IDLE=s        after the solver is compiled, sit idle s seconds, then report memory
#   HTG_SELFTEST_SHOT=path     write a PNG screenshot of the window there

const SELFTEST_FAILED = Ref(false)

const SELFTEST_JS = raw"""
(async () => {
  const out = [];
  const sleep = ms => new Promise(r => setTimeout(r, ms));
  const $ = (s, root = document) => root.querySelector(s);
  const $$ = (s, root = document) => Array.from(root.querySelectorAll(s));
  const visible = el => !!el && el.offsetParent !== null;
  const check = (name, ok, detail) => out.push({name, ok: !!ok, detail: String(detail)});
  const click = el => { el.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
                        el.dispatchEvent(new MouseEvent('mouseup', {bubbles: true})); el.click(); };
  const waitFor = async (cond, ms = 10000) => {
    const t = performance.now();
    while (performance.now() - t < ms) { if (cond()) return Math.round(performance.now() - t); await sleep(25); }
    return -1;
  };
  try {
    await document.fonts.ready;
    check('page title', document.title === 'HTG equilibrium', document.title);
    const header = $('.q-header');
    const title = $('.q-header .text-h6');
    const status = $('.q-header .htg-status');
    check('fixed header: title left, status line right', header && title && status &&
          title.textContent.trim() === 'HTG equilibrium' && getComputedStyle(header).position === 'fixed' &&
          title.getBoundingClientRect().left < 40 && status.getBoundingClientRect().right > window.innerWidth - 40,
          'title "' + (title && title.textContent.trim()) + '", status "' + (status && status.textContent.trim()) + '"');
    const bg = getComputedStyle(header).backgroundColor;
    check('header colour (NiceGUI primary #5898d4)', bg === 'rgb(88, 152, 212)', bg);

    // tabs: names, icons as glyphs (the bundled font), not as words
    const tabs = $$('.q-tab');
    const labels = tabs.map(t => ($('.q-tab__label', t) || {}).innerText);
    const icons = tabs.map(t => ($('.q-icon', t) || {}).textContent);
    check('five tabs with their names', labels.join(',') === 'CASE,QUEUE,RESULTS,COMPARE,SPECIES', labels.join(', '));
    check('tab icons', icons.join(',') === 'edit_note,queue,insights,compare_arrows,science', icons.join(', '));
    const widths = tabs.map(t => Math.round($('.q-icon', t).getBoundingClientRect().width));
    check('icons render as glyphs (Material Icons font loaded)',
          document.fonts.check('24px "Material Icons"') && widths.every(w => w > 0 && w < 30), 'icon widths ' + widths.join(', ') + ' px');

    // tab switching
    let switched = [];
    for (const name of ['Queue', 'Results', 'Compare', 'Species', 'Case']) {
      click(tabs.find(t => t.textContent.includes(name)));
      const w = await waitFor(() => { const p = $$('.q-tab-panel').filter(visible); return p.length === 1 && p[0].textContent.includes('The ' + name + ' tab'); }, 5000);
      const active = ($('.q-tab--active .q-tab__label') || {}).innerText;
      switched.push(name + (w >= 0 && active === name.toUpperCase() ? ' ok' : ' FAILED'));
    }
    check('tabs switch (panel and active tab follow)', switched.every(s => s.endsWith(' ok')), switched.join(', '));
    const panel = $$('.q-tab-panel').filter(visible)[0];
    check('placeholder text', panel && panel.textContent.trim().startsWith('The Case tab is not in the Julia version yet.'),
          panel && panel.textContent.trim());
    const panels = $('.htg-panels');
    check('panel container: max width 1536 px, centred', getComputedStyle(panels).maxWidth === '1536px', getComputedStyle(panels).maxWidth);

    // nothing from the internet
    const res = performance.getEntriesByType('resource').map(r => r.name);
    const remote = res.filter(u => !u.startsWith(location.origin));
    check('page resources all from this server', remote.length === 0, res.length + ' resources' + (remote.length ? '; remote: ' + remote.join(' ') : ''));
  } catch (e) {
    out.push({name: 'self-test script', ok: false, detail: String(e) + ' ' + e.stack});
  }
  return out;
})()
"""

"""The status line text shown in the window now."""
window_status(win) = run(win, "(document.querySelector('.q-header .htg-status') || {}).textContent || ''")

function wait_until(f; timeout=60, step=0.1)
    t = time()
    while time() - t < timeout
        f() && return true
        sleep(step)
    end
    return false
end

function report!(results, name, ok, detail)
    push!(results, Dict("name" => name, "ok" => ok, "detail" => string(detail)))
end

"""Drive the page and check it (on_usable hook of run_window); closes the window at the end."""
function selftest(app, win)
    results = Any[]
    try
        _selftest(app, win, results)
    catch exc
        report!(results, "self-test", false, sprint(showerror, exc, catch_backtrace()))
    end
    npass = 0
    for r in results
        npass += r["ok"]
        println("[htg]   ", r["ok"] ? "PASS" : "FAIL", "  ", rpad(r["name"], 58), r["detail"])
    end
    SELFTEST_FAILED[] = npass < length(results)
    stamp("self-test: $(npass)/$(length(results)) passed")
    stamp("closing the window as a user would (BrowserWindow.close)")
    ElectronAPI.close(win)
end

function _selftest(app, win, results)
    state = STATE[]
    run(app, "electron.BrowserWindow.fromId($(win.id)).webContents.focus()")
    first_status = window_status(win)
    stamp("self-test: driving the page (status line: \"$first_status\")")
    append!(results, run(win, SELFTEST_JS))

    # settings: read at start-up, written when changed
    expect = get(ENV, "HTG_SELFTEST_EXPECT", "")
    isempty(expect) || report!(results, "setting read at start-up (pressure_unit)",
                               state.settings["pressure_unit"] == expect,
                               "\"$(state.settings["pressure_unit"])\", expected \"$expect\" from $(state.settings.path)")
    set = get(ENV, "HTG_SELFTEST_SET", "")
    if !isempty(set)
        state.settings["pressure_unit"] = set
        report!(results, "setting written (pressure_unit)", HTGEquilibrium.load_settings(dirname(state.settings.path))["pressure_unit"] == set,
                "\"$set\" in $(state.settings.path)")
    end

    if get(ENV, "HTG_SELFTEST_EARLY", "") == "1"
        # close the window during the loading or the warm-up: the process must end at once
        report!(results, "a worker is still loading or compiling when the window closes", occursin(r"Loading|starting|warming", status_text(state)),
                "\"$(status_text(state))\"")
        print_processes(app)
        return results
    end

    # the model loads and the pool starts behind the window
    loaded = wait_until(() -> ready(state) || !isempty(state.problem); timeout=300)
    stamp("model loaded, pool started (status line: \"$(window_status(win))\")")
    report!(results, "model loaded and pool started in the background", loaded && isempty(state.problem),
            isempty(state.problem) ? "$(length(state.feed_species)) feed species, $(state.pool.target) workers" : state.problem)
    report!(results, "status line while loading", first_status == "Loading the species data" || startswith(first_status, "Workers: "),
            "\"$first_status\"")

    # the status line follows the pool (server push every second)
    if loaded && isempty(state.problem)
        n = state.pool.target
        wait_until(() -> startswith(window_status(win), "Workers: "); timeout=5)
        before = window_status(win)
        set_workers(state.pool, n + 1)
        t = time()
        changed = wait_until(() -> window_status(win) != before; timeout=5)
        after = window_status(win)
        set_workers(state.pool, n)
        report!(results, "status line updates (set_workers $n -> $(n + 1))", changed,
                "\"$before\" -> \"$after\" after $(round(Int, 1000 * (time() - t))) ms")
        warm = wait_until(() -> !any(w -> w["state"] in ("warming", "starting"), status(state.pool)["workers"]); timeout=300)
        stamp("solver compiled: warm-up done")
        # the status line follows within a second
        shown = wait_until(() -> !occursin(r"warming|starting", window_status(win)); timeout=5)
        report!(results, "warm-up finished (status line follows)", warm && shown, "\"$(window_status(win))\"")
    end

    blocked = blocked_requests(app)
    report!(results, "no request left this computer (Electron filter)", isempty(blocked),
            isempty(blocked) ? "none blocked" : join(blocked, " "))
    report!(results, "no .theme file in the working directory", !ispath(joinpath(pwd(), ".theme")), pwd())

    shot = get(ENV, "HTG_SELFTEST_SHOT", "")
    isempty(shot) || screenshot(app, win, shot)

    idle = something(tryparse(Float64, get(ENV, "HTG_SELFTEST_IDLE", "0")), 0.0)
    if idle > 0
        sleep(idle)
        memory_report(app)
    end
    print_processes(app)
    return results
end

print_processes(app) =
    Sys.islinux() && println("[htg] processes: julia $(getpid()); electron $(join(process_tree(getpid(app.proc)), " "))")

"""Write a PNG of the window to `path`."""
function screenshot(app, win, path)
    js = """electron.BrowserWindow.fromId($(win.id)).webContents.capturePage()
              .then(img => { require('fs').writeFileSync($(JSON.json(path)), img.toPNG()); return true })"""
    electron_await(app, js)
    stamp("screenshot written to $path")
end

# ---------------------------------------------------------------- memory (Linux, M0 spike)
function children_of(pid)
    kids = Int[]
    for d in readdir("/proc")
        all(isdigit, d) || continue
        stat = try read("/proc/$d/stat", String) catch; continue end
        fields = split(stat[findlast(')', stat)+2:end])
        parse(Int, fields[2]) == pid && push!(kids, parse(Int, d))
    end
    return kids
end

process_tree(pid) = vcat([pid], [p for k in children_of(pid) for p in process_tree(k)])

function mem_kb(pid, key)
    file = "/proc/$pid/smaps_rollup"
    isfile(file) || return 0
    for line in eachline(file)
        startswith(line, key * ":") && return parse(Int, split(line)[2])
    end
    return 0
end

function memory_report(app)
    Sys.islinux() || return nothing
    jl = getpid()
    el = process_tree(getpid(app.proc))
    jl_rss, jl_pss = mem_kb(jl, "Rss"), mem_kb(jl, "Pss")
    el_rss, el_pss = sum(mem_kb(p, "Rss") for p in el), sum(mem_kb(p, "Pss") for p in el)
    mb(x) = round(Int, x / 1024)
    println("[htg] memory: julia RSS $(mb(jl_rss)) MB (PSS $(mb(jl_pss)) MB); " *
            "electron $(length(el)) processes RSS $(mb(el_rss)) MB (PSS $(mb(el_pss)) MB); " *
            "total PSS $(mb(jl_pss + el_pss)) MB")
    flush(stdout)
end
