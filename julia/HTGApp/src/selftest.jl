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
#
# The Case tab: the WP4 feed typed into the form (numbers, the ash error), an
# example case opened, changed and saved with Save as (the native dialog is
# answered by the test), the page reloaded, and a case added to the queue.

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

    // tab switching (Case: the case editor, or its "Starting" line; the others: placeholders)
    let switched = [];
    const shows = (p, name) => name === 'Case' ? !!$('[data-mark="case-tab"], .htg-case-starting', p)
                                               : p.textContent.includes('The ' + name + ' tab is not in the Julia version yet.');
    for (const name of ['Queue', 'Results', 'Compare', 'Species', 'Case']) {
      click(tabs.find(t => t.textContent.includes(name)));
      const w = await waitFor(() => { const p = $$('.q-tab-panel').filter(visible); return p.length === 1 && shows(p[0], name); }, 5000);
      const active = ($('.q-tab--active .q-tab__label') || {}).innerText;
      switched.push(name + (w >= 0 && active === name.toUpperCase() ? ' ok' : ' FAILED'));
    }
    check('tabs switch (panel and active tab follow)', switched.every(s => s.endsWith(' ok')), switched.join(', '));
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

# ---------------------------------------------------------------- the Case tab
# Each step is an async script; `body` uses the helpers below and returns nothing
# (the checks go to `out`). Fields are typed as a user types: the input element's
# value is set and an input event sent, so Quasar's q-input sends the change.
const CASE_JS_HEAD = raw"""
(async () => {
  const out = [];
  const sleep = ms => new Promise(r => setTimeout(r, ms));
  const $ = (s, root = document) => root.querySelector(s);
  const $$ = (s, root = document) => Array.from(root.querySelectorAll(s));
  const check = (name, ok, detail) => out.push({name, ok: !!ok, detail: String(detail)});
  const click = el => { el.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
                        el.dispatchEvent(new MouseEvent('mouseup', {bubbles: true})); el.click(); };
  const waitFor = async (cond, ms = 10000) => {
    const t = performance.now();
    while (performance.now() - t < ms) { try { if (cond()) return Math.round(performance.now() - t); } catch (e) {} await sleep(25); }
    return -1;
  };
  const M = () => window.GENIEMODEL;
  const input = key => $('[data-key="' + key + '"]');
  const field = key => input(key) && input(key).closest('.q-field');
  const fieldError = key => { const f = field(key); const m = f && $('.q-field__messages', f);
                              return f && f.classList.contains('q-field--error') && m ? m.textContent.trim() : ''; };
  const type = async (key, text) => {
    const el = input(key);
    el.focus(); el.value = text; el.dispatchEvent(new Event('input', {bubbles: true}));
    await waitFor(() => M().case_texts[key] === text, 2000);
  };
  const settle = async () => { const s = JSON.stringify(M().case_view); await sleep(300);
                               await waitFor(() => JSON.stringify(M().case_view) === s, 3000); };
  const text = mark => ($('[data-mark="' + mark + '"]') || {}).textContent?.trim() ?? '';
  const button = mark => $('[data-mark="' + mark + '"]');
  const disabled = mark => { const b = button(mark); return !!b && (b.disabled || b.classList.contains('disabled')); };
  try {
"""
const CASE_JS_TAIL = raw"""
  } catch (e) {
    out.push({name: 'Case tab script', ok: false, detail: String(e) + ' ' + e.stack});
  }
  return out;
})()
"""
case_js(body) = CASE_JS_HEAD * body * CASE_JS_TAIL

# The WP4 feed (spec WP4 acceptance): numbers on every keystroke; ash 2 kg is more than the solids.
const CASE_FEED_JS = case_js(raw"""
    click($$('.q-tab').find(t => t.textContent.includes('Case')));
    const ready = await waitFor(() => M().case_view.ready && input('feed.total'), 300000);
    check('Case tab: the form is shown', ready >= 0, 'after ' + ready + ' ms');
    M().caseOp('new');
    await waitFor(() => M().case_texts['name.name'] === 'new_case');
    for (const [key, value] of [['feed.total', '10'], ['feed.solids', '15'], ['feed.ash', '8.45'], ['feed.element:C', '30.4'],
                                ['feed.element:H', '4.41'], ['feed.element:O', '53.58'], ['feed.element:N', '2.66'],
                                ['feed.element:S', '0.51']]) await type(key, value);
    // The new case's element percentages are WP4's already: wait for the typed total, solids and
    // ash to reach the Feed summary (the first edits are slow while the solver's warm-up compiles).
    const summaryCell = name => { const c = $$('div', $('[data-mark="case-feed-summary"]')).map(d => d.textContent.trim());
                                  const i = c.indexOf(name); return i < 0 ? '' : c[i + 1]; };
    await waitFor(() => summaryCell('Water') === '8.500' && summaryCell('Ash') === '0.127', 120000);
    await settle();
    const used = ['C', 'H', 'O', 'N', 'S'].map(e => $('[data-used="' + e + '"]').textContent.trim());
    check('Case tab: WP4 feed, Used column', used.join(' ') === '33.20 4.82 58.52 2.91 0.56', used.join(' '));
    check('Case tab: WP4 feed, entered sum and notice', text('case-entered-sum') === '91.56' &&
          text('case-scale-notice') === 'The element percentages add up to 91.56 %. They have been scaled to 100 % (the Used column).',
          text('case-entered-sum') + ' / ' + text('case-scale-notice'));
    const cells = $$('div', $('[data-mark="case-feed-summary"]')).map(d => d.textContent.trim());
    const row = name => { const i = cells.indexOf(name); return i < 0 ? [] : cells.slice(i + 1, i + 3); };
    const rows = {Water: row('Water'), Ash: row('Ash'), 'Ash-free solids': row('Ash-free solids'), C: row('C')};
    check('Case tab: WP4 feed summary (water 8.500 kg 471.8 mol, ash 0.127, ash-free 1.373, C 0.456 kg 37.96 mol)',
          rows.Water.join(' ') === '8.500 471.8' && rows.Ash[0] === '0.127' && rows['Ash-free solids'][0] === '1.373' &&
          rows.C.join(' ') === '0.456 37.96', JSON.stringify(rows));
    check('Case tab: WP4 feed can be queued', !disabled('case-queue') && text('case-status').startsWith('Ready to run'), text('case-status'));
    // ash in kg: the value is converted, then 2 kg is more than the 1.5 kg of solids
    const kg = $$('button', input('feed.ash').closest('.q-field').parentElement).find(b => b.textContent.trim() === 'kg');
    click(kg);
    const converted = await waitFor(() => M().case_texts['feed.ash'] === '0.12675');
    check('Case tab: ash switched to kg (8.45 % of 1.5 kg = 0.12675 kg)', converted >= 0, M().case_texts['feed.ash']);
    await type('feed.ash', '2');
    await waitFor(() => fieldError('feed.ash') !== '');
    check('Case tab: ash 2 kg: message under the field', fieldError('feed.ash') === "Ash can't be more than the solids.", fieldError('feed.ash'));
    check('Case tab: ash 2 kg: Add to queue, Save and Save as disabled, status names the problem',
          disabled('case-queue') && disabled('case-save') && disabled('case-save-as') &&
          text('case-status') === "✗ Ash can't be more than the solids.", text('case-status'));
    await type('feed.ash', '0.1');
    await waitFor(() => !disabled('case-queue'));
    check('Case tab: ash 0.1 kg: allowed again', !disabled('case-queue') && fieldError('feed.ash') === '', text('case-status'));
""")

# Open an example case, change a value, Save as (the native dialog answered by the test).
case_open_js(path) = case_js("""
    M().caseOp('open', $(JSON.json(path)));
    const opened = await waitFor(() => M().case_texts['name.name'] === 'yakaboylu2013_fig7_8' &&
                                       text('case-file-status').includes('yakaboylu2013_fig7_8.toml'), 60000);
    check('Case tab: example case opened', opened >= 0 && input('solver.starts').value === '1', text('case-file-status'));
    await type('solver.starts', '3');
    await waitFor(() => text('case-file-status').endsWith('unsaved changes'));
    check('Case tab: a change is an unsaved change', text('case-file-status').endsWith('· unsaved changes'), text('case-file-status'));
    click(button('case-save-as'));
    await waitFor(() => !text('case-file-status').endsWith('unsaved changes'), 10000);
    check('Case tab: saved with Save as', !text('case-file-status').endsWith('unsaved changes'), text('case-file-status'));
""")

# After the page was reloaded: the saved case is still the one in the form.
const CASE_RELOADED_JS = case_js(raw"""
    const ready = await waitFor(() => M() && M().case_view.ready && input('solver.starts'), 60000);
    check('Case tab after reloading the page: the saved case', ready >= 0 && input('solver.starts').value === '3' &&
          M().case_texts['name.name'] === 'yakaboylu2013_fig7_8' && !text('case-file-status').endsWith('unsaved changes'),
          'Initial guesses ' + (input('solver.starts') || {}).value + ', ' + text('case-file-status'));
""")

# Add to queue: the Queue tab opens.
const CASE_QUEUE_JS = case_js(raw"""
    M().caseOp('new');
    await waitFor(() => M().case_texts['name.name'] === 'new_case' && !disabled('case-queue'));
    click(button('case-queue'));
    const w = await waitFor(() => ($('.q-tab--active .q-tab__label') || {}).innerText === 'QUEUE', 10000);
    check('Case tab: Add to queue opens the Queue tab', w >= 0, ($('.q-tab--active .q-tab__label') || {}).innerText);
""")

"""The Case tab in the window: the WP4 feed, open, change, Save as, reload; Add to queue when the pool runs."""
function case_selftest(app, win, results; queue)
    state = STATE[]
    append!(results, run(win, CASE_FEED_JS))
    # Save as: the native dialog is answered with a file in a new temporary folder
    dir = mktempdir()
    target = joinpath(dir, "saved_case.toml")
    asked = Ref{Any}(nothing)
    NATIVE_SAVE[] = opts -> (asked[] = opts; target)
    example = joinpath(pkgdir(HTGEquilibrium), "..", "..", "cases", "yakaboylu2013_fig7_8.toml")
    append!(results, run(win, case_open_js(normpath(example))))
    opts = asked[]
    report!(results, "Case tab: Save as asks with the native dialog",
            opts !== nothing && opts["title"] == "Save case file as" && endswith(opts["defaultPath"], "yakaboylu2013_fig7_8.toml"),
            opts === nothing ? "not asked" : "\"$(opts["title"])\", $(opts["defaultPath"])")
    saved = isfile(target) ? read(target, String) : ""
    report!(results, "Case tab: the file is written (format 2, starts = 3) and loads",
            occursin("format = 2", saved) && occursin("starts = 3", saved) && load_case(target).starts == 3, target)
    run(win, "location.reload(); true")
    sleep(1)
    wait_until(timeout=60) do
        try
            run(win, "!!(window.GENIEMODEL && window.GENIEMODEL.case_view && window.GENIEMODEL.case_view.ready)") === true
        catch
            false   # the page is loading
        end
    end
    append!(results, run(win, CASE_RELOADED_JS))
    if queue
        before = isdir(state.results_dir) ? readdir(state.results_dir) : String[]
        append!(results, run(win, CASE_QUEUE_JS))
        new_runs() = filter(d -> startswith(d, "new_case-"), setdiff(isdir(state.results_dir) ? readdir(state.results_dir) : String[], before))
        finished = wait_until(() -> !isempty(new_runs()) && all(w -> w["state"] == "idle", status(state.pool)["workers"]); timeout=120)
        report!(results, "Case tab: the queued case runs (5 points, saved below the results folder)", finished,
                join(new_runs(), ", "))
    end
    rm(dir; recursive=true, force=true)
    return results
end

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

    case_selftest(app, win, results; queue=loaded && isempty(state.problem))

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
