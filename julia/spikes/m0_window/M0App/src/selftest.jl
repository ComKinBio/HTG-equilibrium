# Drives the page in the Electron window with JavaScript (clicks, typing, events) and
# checks what the DOM shows afterwards. One line per component: PASS/FAIL and details.
# Native OS dialogs are not driven here (they block on the user); see README.

const SELFTEST_JS = raw"""
(async () => {
  const out = [];
  const sleep = ms => new Promise(r => setTimeout(r, ms));
  const $ = (s, root = document) => root.querySelector(s);
  const $$ = (s, root = document) => Array.from(root.querySelectorAll(s));
  const visible = el => !!el && el.offsetParent !== null;
  const byText = (sel, text) => $$(sel).find(e => e.textContent.trim().includes(text));
  const check = (name, ok, detail) => out.push({name, ok: !!ok, detail: String(detail)});
  const click = el => { el.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
                        el.dispatchEvent(new MouseEvent('mouseup', {bubbles: true})); el.click(); };
  const type = async (input, text) => {
    input.focus(); click(input);
    input.value = text; input.dispatchEvent(new Event('input', {bubbles: true}));
    await sleep(400);
  };
  const tab = async name => { click(byText('.q-tab', name)); await sleep(400); };
  const notes = () => $$('.q-notification').map(n => n.textContent.trim());
  // wait until cond() is true; returns the time taken in ms, or -1 after `ms`
  const waitFor = async (cond, ms = 30000) => {
    const t = performance.now();
    while (performance.now() - t < ms) { if (cond()) return Math.round(performance.now() - t); await sleep(25); }
    return -1;
  };
  try {
    await document.fonts.ready;
    const ic = $('.q-tab .q-icon');
    check('icons (Material Icons font)', ic && ic.getBoundingClientRect().width < 30,
          'tab icon width ' + (ic && ic.getBoundingClientRect().width) + ' px');

    // status line pushed by a server timer
    const st = () => byText('.q-header .text-body2', 'Workers').textContent;
    const s1 = st(); await sleep(1600); const s2 = st();
    check('server push (timer -> status line)', s1 !== s2, s1 + ' -> ' + s2);

    // select with search
    const sel = $$('.q-select')[0];
    const inp = $('input', sel);
    await type(inp, 'alkane C1');
    const opts = $$('.q-menu .q-item').map(e => e.textContent.trim());
    check('select with search: filters while typing', opts.length > 0 && opts.every(o => o.includes('alkane C1')),
          opts.length + ' options for "alkane C1", first: ' + opts[0]);
    click($$('.q-menu .q-item')[0]);
    const wp = await waitFor(() => inp.value.includes('alkane C1'), 10000);
    check('select with search: picks an option', inp.value.includes('alkane C1'), 'input shows "' + inp.value + '" after ' + wp + ' ms');
    inp.blur(); document.body.click(); await sleep(200);

    // multi-select with chips
    const msel = $$('.q-select')[1];
    const minp = $('input', msel);
    await type(minp, 'fig7');
    await waitFor(() => $$('.q-menu .q-item').length > 0, 10000);
    let mopts = $$('.q-menu .q-item');
    const nfig7 = mopts.length;
    if (mopts.length > 0) { click(mopts[0]); await waitFor(() => $$('.q-chip', msel).length === 1, 10000); }
    mopts = $$('.q-menu .q-item'); if (mopts.length > 1) { click(mopts[1]); await waitFor(() => $$('.q-chip', msel).length === 2, 10000); }
    const chipsIn = $$('.q-chip', msel).map(c => c.textContent.trim());
    check('multi-select with chips + search', chipsIn.length === 2, nfig7 + ' matches; chips: ' + chipsIn.join(' | '));
    minp.blur(); document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', keyCode: 27, bubbles: true}));
    document.body.click(); await sleep(300);

    // number field with unit switch
    const tIn = $$('input[type=number]')[0];
    const t0 = tIn.value;
    click(byText('.q-btn-toggle .q-btn', 'K'));
    const w1 = await waitFor(() => tIn.value == '623.15');
    const t1 = tIn.value;
    click(byText('.q-btn-toggle .q-btn', '°C'));
    const w2 = await waitFor(() => tIn.value == '350');
    const t2 = tIn.value;
    check('number field + unit switch (°C/K)', t0 == '350' && t1 == '623.15' && t2 == '350',
          t0 + ' °C -> ' + t1 + ' K -> ' + t2 + ' °C; server round trip ' + w1 + ' ms (first), ' + w2 + ' ms (second)');
    tIn.focus(); tIn.value = ''; tIn.dispatchEvent(new Event('input', {bubbles: true})); tIn.blur(); await sleep(300);
    const rule = $$('.q-field--error').length;
    tIn.focus(); tIn.value = '350'; tIn.dispatchEvent(new Event('input', {bubbles: true})); tIn.blur(); await sleep(300);
    check('number field validation rule', rule > 0, rule + ' field(s) in error state when empty');

    // tooltip (live text from the server)
    const gIn = $$('input[type=number]')[1];
    gIn.focus(); gIn.value = '7'; gIn.dispatchEvent(new Event('input', {bubbles: true}));
    await waitFor(() => Genie.findApp().starts_tooltip.startsWith('7'));
    const hover = () => { gIn.closest('.q-field').dispatchEvent(new MouseEvent('mouseenter', {bubbles: false}));
                          gIn.closest('.q-field').dispatchEvent(new MouseEvent('mouseover', {bubbles: true})); };
    for (let i = 0; i < 5 && !$('.q-tooltip'); i++) { hover(); await waitFor(() => $('.q-tooltip'), 1500); }
    await sleep(200);
    const tips = $$('.q-tooltip').map(t => t.textContent.trim());
    check('tooltip (with live text)', tips.some(t => t.includes('7 initial guesses')), 'tooltips: ' + tips.join(' | '));

    // option group, radio, checkbox, chips
    click(byText('.q-radio', 'SLSQP (q-radio)')); await sleep(400);
    const og = $$('.q-option-group .q-radio__inner--truthy').length;
    check('radio + option group (shared value)', og === 1 && byText('.q-option-group .q-radio', 'SLSQP').querySelector('.q-radio__inner--truthy'),
          'option group follows the radio');
    click(byText('.q-checkbox', 'Continuation')); await sleep(300);
    check('checkbox', $('.q-checkbox .q-checkbox__inner--falsy'), 'unticked after click');
    const aq = byText('.q-chip', 'Aqueous');
    click(aq); await sleep(400);
    check('selectable chips', aq.classList.contains('q-chip--selected'), aq.className);
    check('badge', byText('.q-badge', '3 running'), '');

    // expansion item
    click(byText('.q-expansion-item .q-item', 'Advanced settings')); await sleep(600);
    check('expansion item', visible($('.q-expansion-item textarea')), 'textarea visible after opening');
    check('CSS grid (ui.grid)', getComputedStyle($('.m0-grid')).display === 'grid', '');

    // notifications
    click(byText('.q-btn', 'Notify'));
    const wn = await waitFor(() => notes().some(n => n.includes('Saved case file')));
    check('notification from the server', wn >= 0, 'after ' + wn + ' ms: ' + notes().join(' | '));
    click(byText('.q-btn', 'Check (spinner'));
    await waitFor(() => $$('.q-notification .q-spinner').length > 0);
    const spin = $$('.q-notification .q-spinner').length;
    await waitFor(() => notes().some(n => n.includes('All points checked')));
    await waitFor(() => $$('.q-notification .q-spinner').length === 0, 5000);
    check('persistent notification with spinner, dismissed by server', spin > 0 && $$('.q-notification .q-spinner').length === 0 && notes().some(n => n.includes('All points checked')),
          'spinner seen: ' + spin + '; after: ' + notes().join(' | '));

    // dialog
    click(byText('.q-btn', 'Open dialog')); await sleep(600);
    const dlg = byText('.q-dialog', 'Remove species?');
    click(byText('.q-dialog .q-btn', 'Remove'));
    await waitFor(() => notes().some(n => n.includes('Dialog confirmed'))); await sleep(500);
    check('dialog (open, button handled on server, closed by server)', dlg && !byText('.q-dialog', 'Remove species?') && notes().some(n => n.includes('Dialog confirmed')),
          'opened: ' + !!dlg);

    // fixed header
    window.scrollTo(0, document.body.scrollHeight); await sleep(300);
    const top = $('.q-header').getBoundingClientRect().top;
    check('fixed header with tabs while scrolling', window.scrollY > 200 && top === 0, 'scrollY ' + window.scrollY + ', header top ' + top);
    window.scrollTo(0, 0);

    // queue tab: linear progress, list, scroll area
    await tab('Queue');
    check('tabs + tab panels', visible($('.q-linear-progress')), 'Queue panel visible');
    check('linear progress, list, badges', $$('.q-list .q-item').length === 2 && byText('.q-badge', 'finished'), '');
    check('scroll area with log text', byText('.q-scrollarea pre', 'point 1/31'), '');

    // results tab: plot, table, image, unit switch redraws plot
    await tab('Results'); await sleep(800);
    const traces = $$('.js-plotly-plot .scatterlayer .trace').length;
    const xt = () => ($('.js-plotly-plot .xtitle') || {}).textContent;
    const x0 = xt();
    click(byText('.q-tab-panel .q-btn-toggle .q-btn', 'K'));
    await waitFor(() => (xt() || '').includes('(K)'));
    const x1 = xt();
    click(byText('.q-tab-panel .q-btn-toggle .q-btn', '°C'));
    await waitFor(() => (xt() || '').includes('°C'));
    check('StipplePlotly plot (3 traces, redrawn on unit switch)', traces === 3 && x0.includes('°C') && x1.includes('K'), traces + ' traces; x title "' + x0 + '" -> "' + x1 + '"');
    check('table', $$('.q-table tbody tr').length === 12, $$('.q-table tbody tr').length + ' rows, columns: ' + $$('.q-table thead th').map(t => t.textContent.trim()).join(', '));
    const img = $('.q-img img');
    check('image (q-img, data URL)', img && img.complete, img ? img.src.slice(0, 30) : 'none');

    // files tab: Quasar uploader driven through its API (the OS picker itself is Chromium's)
    await tab('Files');
    const vm = Genie.findApp();
    const up = vm.$refs.uploader;
    up.addFiles([new File(['name,formula\nmethane,CH4\n'], 'species.csv', {type: 'text/csv'})]);
    const wu = await waitFor(() => (Genie.findApp().upload_info || '').includes('species.csv'));
    const info = byText('.q-tab-panel div', 'uploaded');
    check('uploader (file content to server)', info && info.textContent.includes('species.csv: 25 bytes'),
          'after ' + wu + ' ms: ' + (info ? info.textContent.trim() : 'no upload info'));

    await tab('Case');
  } catch (e) {
    out.push({name: 'selftest script', ok: false, detail: String(e) + ' ' + e.stack});
  }
  return out;
})()
"""

# Real (OS-level) mouse and keyboard events through Electron, as a user would type:
# click into the species select and type "phen"; the filtered menu must show phenol.
function real_keyboard_check(app, win)
    wc = "electron.BrowserWindow.fromId($(win.id)).webContents"
    run(app, "$wc.focus()")
    r = run(win, "(() => { const b = document.querySelectorAll('.q-select')[0].querySelector('input').getBoundingClientRect();" *
                 " return [Math.round(b.x + 20), Math.round(b.y + b.height / 2)] })()")
    x, y = r
    run(app, "($wc.sendInputEvent({type: 'mouseDown', x: $x, y: $y, button: 'left', clickCount: 1}), " *
             "$wc.sendInputEvent({type: 'mouseUp', x: $x, y: $y, button: 'left', clickCount: 1}), null)")
    sleep(0.3)
    for ch in "phen"
        run(app, "$wc.sendInputEvent({type: 'char', keyCode: '$ch'})")
    end
    sleep(0.8)
    opts = run(win, "Array.from(document.querySelectorAll('.q-menu .q-item')).map(e => e.textContent.trim())")
    run(app, "($wc.sendInputEvent({type: 'keyDown', keyCode: 'Escape'}), $wc.sendInputEvent({type: 'keyUp', keyCode: 'Escape'}), null)")
    run(win, "(() => { const i = document.querySelectorAll('.q-select')[0].querySelector('input'); i.blur(); return 1 })()")
    sleep(0.3)
    Dict("name" => "select with search: real keyboard input", "ok" => opts == ["phenol"],
         "detail" => "typed \"phen\" -> options $(opts)")
end

function selftest(app, win)
    stamp("self-test: driving the page")
    run(app, "electron.BrowserWindow.fromId($(win.id)).webContents.focus()")  # DOM focus events need it
    results = Any[run(win, SELFTEST_JS)...]
    sleep(0.5)
    push!(results, real_keyboard_check(app, win))
    npass = 0
    for r in results
        npass += r["ok"]
        println("[M0]   ", r["ok"] ? "PASS" : "FAIL", "  ", rpad(r["name"], 58), r["detail"])
    end
    stamp("self-test: $(npass)/$(length(results)) passed")
    results
end
