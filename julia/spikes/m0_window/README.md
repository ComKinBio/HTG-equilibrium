# M0 spike: Genie/Stipple in an Electron.jl window

Throw-away app for milestone M0 of the Phase 2 plan (`docs/spec/julia_migration_plan.md`,
"Window and toolkit test"). It is not part of the app and can be deleted once the
decision is made.

**Result (Linux): go.** With a sysimage the window is usable 2.2–2.4 s after launch
(criterion: under ~20 s), all 26 component checks pass, and closing the window ends
every process. Windows and macOS were out of scope for this run and are still open
(see "Recommendation").

## What it is

| File | Purpose |
|---|---|
| `M0App/src/ui.jl` | The page: Stipple model, handlers and view with the Quasar components Phase 1 uses |
| `M0App/src/shell.jl` | Genie server on a free local port, Electron.jl window, native dialogs, clean shut-down, memory report |
| `M0App/src/selftest.jl` | Drives the page in the window (clicks, typing, real keyboard input) and checks each component |
| `run.jl` | Launcher |
| `precompile.jl`, `build_sysimage.jl` | Sysimage with PackageCompiler; the workload is a full self-test run |
| `measure.sh` | Start-up time and idle memory, N runs without and with the sysimage |

`M0App` is a package (developed into this environment) so that the sysimage contains
the page code the way it would contain the real app. Versions: Julia 1.12.7, Genie 6.4.1,
Stipple 1.3.1, StippleUI 1.1.0 (Quasar 2, Vue 3), StipplePlotly 1.1.2 (plotly.js),
Electron.jl 7.0.0 (Electron 34.2.0), PackageCompiler 2.4.5. All pinned in `Manifest.toml`.

## How to run

```sh
cd julia/spikes/m0_window
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # first time; downloads Electron (263 MB)
julia --project=. run.jl                               # opens the window; close it to quit

julia --project=. build_sysimage.jl                    # ~10 min, opens the window once
julia --project=. -J m0_sysimage.so run.jl             # fast start

M0_SELFTEST=1 M0_AUTOCLOSE=1 julia --project=. -J m0_sysimage.so run.jl   # check every component
bash measure.sh 3                                      # timing and memory
```

Other switches are listed at the top of `run.jl` (`M0_IDLE`, `M0_SHOT`, `M0_GTK_DIALOGS`,
`M0_JSDIR`). Needs a desktop session (Wayland or X11). The sysimage (823 MB) is not
committed.

## Measurements

Machine: Intel i7-1370P (20 threads), 29 GB RAM, Ubuntu 26.04, Linux 7.0, GNOME on
Wayland (Electron runs through XWayland). Another Julia job (the M1 test) was running
at the same time, load average 1.8–3.4 during the timed runs, so the no-sysimage
numbers are somewhat pessimistic.

Time from launch (shell timestamp before `julia` starts) to a **usable** window, meaning
the page is loaded, Vue is mounted, the websocket is connected and a server handler has
run for the client:

| | run 1 | run 2 | run 3 |
|---|---|---|---|
| No sysimage, packages precompiled ("cold") | 41.7 s | 46.9 s | 45.3 s |
| Sysimage | 2.36 s | 2.39 s | 2.24 s |

Without the sysimage, loading the packages takes ~3.6 s and Electron starts in 0.3 s; the
rest (33–38 s) is Julia compiling Genie/Stipple's page rendering for the first request.
Every first interaction is also slow without the sysimage. With the sysimage the page is
shown 1.6–2.4 s after launch and nothing is compiled later.

Idle memory, 10 s after start-up (Linux `smaps_rollup`; RSS counts shared pages in every
process, PSS splits them, so PSS is the fair total):

| | Julia RSS | Electron (6 processes) RSS | Julia + Electron PSS |
|---|---|---|---|
| No sysimage | 1.32 GB | 0.66 GB | 1.67 GB |
| Sysimage | 0.82 GB | 0.66 GB | 1.16 GB |

The Julia figure has no science packages yet (Clapeyron, Ipopt); expect more once they
are in the sysimage.

Server round trip for a field change (unit switch, notification) is ~330 ms, of which
300 ms is Stipple's default debounce on client-to-server field updates (configurable
per field with `@debounce`); an upload round trip was 25 ms.

Building the sysimage took 9.5 min with a peak of 15.7 GB RAM (PackageCompiler) and
gives an 823 MB file. The workload opens the Electron window, so CI needs a display
(xvfb on Linux).

## Component checklist

Phase 1 uses (counted in `src/htg/gui/`) labels, rows/columns/grids, cards, buttons,
selects (with search, multiple, chips), inputs, number fields, textareas, toggles,
radios, checkboxes, tooltips, dialogs, tables, notifications (also persistent with a
spinner), tabs in a fixed header, expansion items, selectable chips, badges, icons,
linear progress, list/items, a log, an image, upload, download, native file dialogs,
timers and Plotly plots. Each was built and checked in the Electron window by the
self-test (26/26 pass with and without the sysimage, also with no network).

| Phase 1 (NiceGUI) | Stipple / Quasar | Status |
|---|---|---|
| `ui.header` with `ui.tabs`, `ui.tab_panels` | `quasar(:layout)`, `quasar(:header)` (raw), `tabgroup`, `tab`, `tabpanels`, `tabpanel` | works; header stays fixed while scrolling |
| `ui.select(with_input=True)` | `select(...; useinput, @on(:filter, "filterSpecies"))` | works, with a ~5-line JS filter method per select (Quasar filters client-side via `@filter`); real keyboard typing checked |
| `ui.select(multiple=True, use-chips)` | `select(...; multiple, usechips, useinput)` | works |
| `ui.number` + unit `ui.toggle` | `numberfield`, `btntoggle`; conversion in `@onchange` | works; validation via `:rules` |
| `ui.input`, `ui.textarea` | `textfield`, `textarea` | works |
| `ui.radio` (is a q-option-group) | `quasar(:option__group)` (raw) or `radio` | works |
| `ui.checkbox` | `checkbox` | works |
| `ui.tooltip` (live text) | `tooltip("{{ field }}")` | works |
| `ui.dialog` | `dialog(:show; persistent)` | works; buttons handled on the server |
| `ui.table` | `table(:t)` with `DataTable(DataFrame)` | works |
| `ui.plotly` | StipplePlotly `plot(:data, layout = :layout)` | works; redraws on data change |
| `ui.expansion` | `expansionitem` | works |
| `ui.chip(selectable=True)` | `chip(; selected = :field)` | works |
| `ui.badge`, `ui.icon` | `badge`, `icon` | works; icon font needs a workaround (below) |
| `ui.notify` | `notify(model, msg, :type)` | works |
| `ui.notification(spinner, timeout=None)` then dismiss | `$q.notify({spinner, timeout: 0})` in `@methods`, started/dismissed with `run(model, "this.startCheck()")` | workaround (small JS) |
| `ui.linear_progress` | `quasar(:linear__progress)` (raw) | works |
| `ui.list`, `ui.item` | `quasar(:list)`, `item`, `itemsection` | works |
| `ui.log` | `quasar(:scroll__area)` + `<pre>` | workaround (no log component) |
| `ui.image` (save-plot preview) | `imageview(src = :field)` with a data URL | works |
| `ui.grid`, `ui.element` swatches | `htmldiv(style = "display:grid; ...")`, `span` | works |
| `ui.timer` (status line) | `@async` loop in Julia setting a field | works (server push) |
| `ui.upload` | `uploader` + `StippleUI.enable_uploads()` (off by default in prod) | works; content arrives as a temp file |
| `ui.download` / `save_bytes` | Electron save dialog, Julia writes the file; or a browser download | both work (see below) |
| native open dialog (pywebview) | Electron `dialog.showOpenDialog` from Julia | works (see below) |
| `ui.refreshable` | Vue `v-if`/`v-for` over reactive fields | not built in the spike; a standard Vue pattern, but every refreshable section needs rewriting |
| `.mark()` test hooks | none; the self-test drives the page with JS in Electron | different test approach for M4 |

Nothing in Phase 1 is missing.

## File pickers

- **Browser picker (QUploader / `<input type=file>`) is not enough on its own.** In
  Electron 32+ the page cannot get a file's path, only its content. That is fine for
  "Import from file" (the content is uploaded to the server), but not for "Save as" or
  anything that needs a path.
- **Electron's native dialogs work from Julia.** `shell.jl` calls
  `electron.dialog.showOpenDialog` / `showSaveDialog` in Electron's main process
  through `run(app, js)`. Electron.jl cannot await a JavaScript Promise, so a small
  helper (`electron_await`) parks the result in a global and polls it. Checked with
  real keystrokes: the open dialog returned the chosen path to the Julia handler, the
  save dialog returned a path and Julia wrote the file there. Filters (e.g.
  "Species files") and titles work.
- On GNOME, Electron uses the XDG desktop portal by default, so the dialog is GNOME's
  own file chooser (checked: a portal request is opened). `M0_GTK_DIALOGS=1` switches to
  Electron's in-process GTK dialog (used for the scripted check, because the portal
  window cannot be driven from a script).
- A browser download (Blob + `<a download>`) also opens Electron's native Save dialog,
  so plain downloads are acceptable too.
- Recommendation: native Electron dialogs for open, save and save-as (as Phase 1 does
  with pywebview); keep the uploader only as a fallback.

## Window and processes

- Electron.jl's `main.js` keeps Electron running after its last window closes, so
  `shell.jl` waits for the window to close, then closes the Electron app, waits for the
  process, stops the Genie server and exits Julia.
- Checked: closing the window (as the window manager does) ends Julia and all 6 Electron
  processes within ~0.3 s (`measure.sh` checks for leftovers after every run: none).
  Killing Julia also ends Electron (its pipe closes), and killing Electron ends Julia.
- The server listens on 127.0.0.1 only, on a free port. The page runs with
  `nodeIntegration = false`, `contextIsolation = true`.

## Problems found and workarounds

1. **Icons need internet by default.** Stipple's default layout loads the Material Icons
   font from fonts.googleapis.com, and its default theme `@import`s the Outfit font from
   there too. Offline, icons show as words (e.g. "edit_note"). Workaround in `ui.jl`:
   `core_theme = false` and serve a local Material Icons woff2 (the spike borrows the copy
   NiceGUI ships). The real app must bundle the font (and use a theme without the
   `@import`). Checked with no network: 26/26 pass.
2. **Name clashes.** `Genie.Renderer.Html` (re-exported by Stipple) has `select`,
   `table`, `textarea`; import the StippleUI ones explicitly
   (`using StippleUI.Selects: select` etc.). A model field must not share a name with a
   function (`plot_layout` broke rendering).
3. **`DEFAULT_LAYOUT(head_content = "...")` splits a String into characters** (it
   `join`s the argument); pass a Vector.
4. **Stipple writes `.theme` into the working directory.** The app should start with a
   working directory it owns (e.g. the user data folder).
5. **Quasar props not known to StippleUI** are passed through literally
   (`fillinput` becomes `fillinput`); use `var"fill-input" = true`.
6. **300 ms debounce** on every client-to-server field change by default; set lower per
   field where it matters.
7. **Uploads are disabled in production mode** until `StippleUI.enable_uploads()`.
8. **Electron.jl ships Electron 34.2.0**, which is past Electron's support window
   (no more security fixes). The app loads only its own local page, so the risk is
   small, but the owner should decide (see below).
9. **Sysimage build is heavy**: 15.7 GB peak RAM, 9.5 min, 823 MB. GitHub's standard
   runners for private repositories have less RAM; check this early in M5.
10. Without a sysimage the app is unusable in practice (40–47 s to start, then a slow
    first click on every control). The sysimage is required, not optional.
11. Electron prints harmless `GetVSyncParametersIfAvailable() failed` lines on XWayland.

## Recommendation

**Go for Genie/Stipple + Electron.jl**, on the Linux evidence: 2.2–2.4 s to a usable
window with a sysimage (criterion ~20 s), every Phase 1 component works (two need a few
lines of JavaScript, two need small layout workarounds), file dialogs are native, and
shut-down is clean. The fallbacks (Bonito.jl in Electron, Gtk4.jl + Gtk4Makie.jl) are
not needed.

The plan's go criterion also asks for Windows and macOS. Before M4 starts, run this
spike there, ideally as the first M5 CI job (build the sysimage and run
`M0_SELFTEST=1 M0_AUTOCLOSE=1` on windows-latest and macos-latest, under xvfb on
ubuntu-latest). Things to watch there: Electron.jl's pipe set-up on Windows, the
macOS app bundle and Gatekeeper, and sysimage build memory on CI runners.

## Decisions for the owner

1. Accept the conditional go (Linux measured; Windows and macOS checked in CI before M4)?
2. Electron version: accept Electron.jl's bundled Electron 34.2.0 for now, or plan in M5
   to supply a current Electron (Electron.jl looks for its artifact; an artifact override
   or a fork of its `Artifacts.toml` would be needed)?
3. File dialogs: use Electron's native dialogs for open/save/save-as (recommended), with
   the GNOME portal dialog on Linux (the default)?
4. Fonts: bundle Material Icons (Apache 2.0) and use a plain system font, so the app
   works offline and makes no requests to Google?
