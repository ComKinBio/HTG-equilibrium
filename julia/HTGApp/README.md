# HTGApp.jl

The Julia version of the HTG equilibrium app (`htg gui`; Phase 2, milestone
M4 of `docs/spec/julia_migration_plan.md`). It is built on the science core
`julia/HTGEquilibrium` and on Genie, Stipple and StippleUI (the page, Vue and
Quasar) and Electron.jl (the window), the stack chosen in the M0 spike
(`julia/spikes/m0_window`).

So far it has the shell of the Python app (`src/htg/gui/app.py`): the window,
the header with the title and the status line, the five tabs, the app state
and the app preferences. The tabs show a placeholder; their screens are ported
in the next steps of M4. Until then the Python app (`uv run htg gui`) is the
one to use.

It is a separate package so that the core (and `htg run`, `check`, ...) does
not load the GUI packages.

## Running it

Once, from the repository root (downloads and precompiles the packages and
Electron, then builds the SUPCRT92 water library of the core if it is not
there yet):

```
julia --project=julia/HTGApp -e 'using Pkg; cd(mktempdir()); Pkg.instantiate()'
julia julia/HTGEquilibrium/deps/build.jl
```

The `cd(mktempdir())` is there because Stipple writes a file `.theme` into
the working folder when it is precompiled (its precompile workload renders a
page with its theme). The app itself never writes it; `.theme` is in the
repository's `.gitignore` in case a precompile runs in the repository.

Then:

```
julia/htg gui                     # Linux, macOS: the app in its own window
julia/htg gui --browser           # in the web browser instead (port 8080)
julia\htg.cmd gui                 # Windows
julia --project=julia/HTGApp -t auto -m HTGApp gui      # without the launcher
```

The options are Python's: `--browser`, `--no-browser` (browser mode: serve
only), `--workers N` (worker threads), `--port`, `--host`, `--results DIR`,
`--cases DIR`, `--user-data DIR`, `--no-user-species`, and `--data DIR`
before `gui`. The launchers start Julia with `--threads=auto` as for `htg run`;
`--workers N` is capped at the threads Julia has, with a note.

The window opens first (a start-up page, then the app page when the server
has compiled it). The species data are loaded and the worker pool is started
in the background; the status line says "Loading the species data" until
then, and the usual "Workers: ..." after. Closing the window stops the
workers and the server and prints "Stopped. Finished runs are saved in
<folder>". Ctrl+C in the terminal does the same.

If the window cannot open (no display, a missing system library, Electron
not installed), one line says why and the app opens in the browser instead.

## Files

| File | Contents | Python counterpart |
|---|---|---|
| `Project.toml`, `Manifest.toml` | The core as a `[sources]` path; Genie 6.4.1, Stipple 1.3.1, StippleUI 1.1.0, Electron 7.0.0 pinned exactly (the M0 versions; added with `PRESERVE_ALL`, so every package the core uses has the core's version) | `pyproject.toml` |
| `Artifacts.toml` | The Material Icons font (the `material-icons` npm package 1.13.14, by sha256), served by the app itself | NiceGUI's bundled font |
| `src/HTGApp.jl` | The module | `gui/__init__.py` |
| `src/state.jl` | `State` (model, pool, settings, results and cases folders, `phase_options`, `feed_species`, `submit`), `use_model!`, `load!` (model and pool in the background), `status_text`, `worker_summary`, `output_folder` | `gui/app.py` (`State`), `gui/common.py` |
| `src/page.jl` | The page: header, status line (pushed every second), tabs and panels; the icon font route; the CSS (NiceGUI's colours, the panel width) | `gui/app.py` (`build_page`) |
| `src/window.jl` | Can the window open? (`window_problem`, `explain_electron_failure`, `window_or_browser`) | `gui/window.py` |
| `src/shell.jl` | `run_app`: server, window, life cycle, native file dialogs (for the tabs), the external-request block | `gui/app.py` (`run_app`) |
| `src/cli.jl` | `main` (`julia -m HTGApp gui ...`; the core's `parse_cli`) | `cli.py` (`cmd_gui`) |
| `src/selftest.jl` | The window self-test (`HTG_SELFTEST=1`) | – |
| `test/runtests.jl` | Unit tests (`Pkg.test`) | `tests/test_gui.py` (shell parts), `tests/test_window.py` |
| `test/dump_state.py`, `test/data/state_python.json` | Python's `phase_options`, `feed_species` and `worker_summary` texts (7.9 kB, committed) | – |
| `test/selftest.jl` | Runs the app three times with the self-test and checks what is left behind | – |

The app preferences (`settings.toml`) are in the core
(`HTGEquilibrium/src/settings.jl`, see the core README): the same file as the
Python app, so both read each other's choices.

## What is ported

- `State`: `model`, `pool`, `settings`, `results_dir`, `cases_dir` (resolved
  as `Path.resolve` does), `phase_options` and `feed_species` (tested equal to
  Python's), `use_model!`, `submit` (into `results/<case>-<yyyymmdd-HHMMSS>`).
  Not yet: `document` (CaseDocument), `comparison_files`, `sources`, `run`,
  `tables`; they come with the tabs that use them.
- The page: title "HTG equilibrium"; the fixed header in the primary colour
  with the title on the left and the status line on the right; the tabs Case
  (edit_note), Queue (queue), Results (insights), Compare (compare_arrows),
  Species (science); the panels below in a container at most 1536 px wide,
  centred (NiceGUI's `max-w-screen-2xl mx-auto`). The colours are NiceGUI's
  (primary #5898d4).
- The status line: `worker_summary(status(pool))`, the same text as Python,
  updated every second.
- `run_app` and `htg gui` with all of Python's options; the window 1400 x 900
  titled "HTG equilibrium"; the browser fallback with one line; "Stopped.
  Finished runs are saved in ..." at the end.

## Differences from Python

- **Window first, data in the background.** Python loaded the species data
  and started the pool, then opened the window. Here the window opens first
  (Julia needs ~10 s to load the packages and ~40 s to compile the page
  without a sysimage, so the user sees something at once), and the data and
  pool are loaded in a background task. Until they are ready the status line
  says "Loading the species data" (Python had no such state); if loading
  fails it says "Could not load the species data: <reason>" and the window
  stays open (Python stopped before opening it). `submit` before that throws
  "The app is still loading the species data. Try again in a moment."
- **A start-up page.** The window shows the header and "Starting the app."
  while the server compiles the real page. Python's window opened when the
  page was ready.
- **Electron instead of pywebview (Qt).** The M0 decision. The window check
  starts a tiny Electron program on Linux (as Python started Qt) and explains
  a missing display or library in one line; the hints name Electron's
  libraries and `Pkg.instantiate` instead of `uv sync`.
- **Threads instead of processes.** `--workers N` means worker threads, as in
  the core's `htg run`: "Starting N worker thread(s) (the first results take
  about a minute while Julia compiles)" instead of "Starting N worker(s); each
  needs about a minute before the first point."
- **Browser mode message.** "HTG equilibrium is running at http://... (Ctrl+C
  stops it)." instead of NiceGUI's welcome message.
- **Placeholders.** Each tab says "The <Tab> tab is not in the Julia version
  yet. For now, use the Python app (uv run htg gui)."
- **Fonts.** The page uses the system's sans-serif font (Roboto if it is
  installed, as in NiceGUI). NiceGUI bundles Roboto; Stipple's default theme
  loads a font from Google, which is turned off here (no internet use). The
  icon font is the npm package's woff2 as a Julia artifact, served by the app.
- **No requests leave the computer.** The window cancels every request that
  is not to the app's own server (fonts, scripts, links to web pages). Python
  did not block anything (NiceGUI served everything itself).
- **Logs.** Genie's and Stipple's start-up messages are hidden (warnings and
  errors only). Warnings use Julia's format (`┌ Warning: ...`).
- **Ctrl+C** in the terminal closes the window and stops the app the same way
  as closing the window (Python: the same, through NiceGUI).
- **The process ends at once.** Julia's normal exit waits until a
  compilation step in progress on a worker thread has finished (20–30 s when
  the window was closed or Ctrl+C pressed during the warm-up). After
  "Stopped. ..." the app therefore ends the process directly (`exit_now`,
  which removes the temporary session folder first; Julia's own exit steps
  are skipped). Python's worker processes ended with the app.
- **Genie sessions.** Genie's session cookie gets a new name in every run and
  no expiry date, and its session files go to a temporary folder of the run:
  the secret token is new in every run, so a 30-day cookie of an earlier run
  could not be read and logged errors on every request.
- **Preferences API.** `settings[key]`, `settings[key] = value`,
  `export_defaults`, `set_export_defaults`; `SettingsKeyError` and
  `SettingsValueError` instead of `KeyError` and `ValueError` (see the core
  README).

## Tests

`julia --project=julia/HTGApp -t 4 -e 'using Pkg; Pkg.test()'`: 70 tests, about 75 s
(2026-11-24).

- `worker_summary` against Python's texts; `output_folder`; path resolution.
- `State`: folders, preferences, `phase_options` and `feed_species` equal to
  Python's, the status line with a pool (a stand-in solver), `submit`.
- Loading in the background: good data, and a broken data folder.
- `gui` options, usage and help; the window check and its messages.
- The page as served by Genie: no reference to the internet, no Google
  font, no Stipple theme files, the font route, and no `.theme` file written
  into the working folder.

**Self-test of the window** (needs a desktop session; Linux for the process
checks):

```
julia julia/HTGApp/test/selftest.jl              # three sessions, ~7 min
julia julia/HTGApp/test/selftest.jl --offline    # the same without network (unshare -rn)
```

It starts the app three times as `htg gui` does, with `HTG_SELFTEST=1`, one
new config folder and an empty working folder each time. In the window it checks: the title, the
fixed header (title left, status line right), the header colour, the tab
names and icons, that the icons are glyphs of the bundled font, that the tabs
switch, the placeholder, the panel width, that every resource came from the
app's server and that no request was blocked; that the data load and the pool
start in the background, and the status line while loading; that the status
line follows the pool (`set_workers`); that the warm-up finishes. The first
session sets the pressure unit to MPa; the second must read it at start-up.
The third closes the window during the warm-up; its process must end within
seconds. After each session: exit code 0, no Julia or Electron process left,
nothing written to the working folder, no error or warning logged, and no
`htg-*` folder left in the system temp folder.

Result (2026-11-24, Linux, Wayland): 18/18, 17/17 and 12/12 checks passed
(the setting was kept), no process, file or log message left, and the third
process ended 0.2 s after its window closed. With `--offline` (no
network at all, loopback only) the same checks pass: the page has 14
resources, all from the app's server, and nothing was blocked.

A manual check without a display (`DISPLAY` and `WAYLAND_DISPLAY` unset)
printed the one line "Cannot open a separate window: no display was found.
Run htg gui from a desktop session. Opening the app in the browser instead
(htg gui --browser does this directly)." and served the page; Ctrl+C
(SIGINT) during the warm-up printed "Stopped. ..." and the process ended
0.2 s later with exit code 0, leaving nothing in the working folder.

## Measurements

Linux, Intel i7-1370P (20 threads), 30 GB, no sysimage, packages
precompiled, two sessions of the self-test (time from launch, before
`julia` starts):

| | Session 1 | Session 2 |
|---|---|---|
| Julia running, HTGApp loaded | 12.6 s | 6.3 s |
| Window open with the start-up page | 20.5 s | 14.1 s |
| Window usable (page loaded, websocket connected) | 62.3 s | 56.1 s |
| Species data loaded, pool started (as seen by the self-test) | 78.8 s | 56.7 s |
| Warm-up done (solver compiled, "Workers: 19 idle") | 179 s | 176 s |

Idle memory after the warm-up (session 2): Julia RSS 2.53 GB (PSS 2.45 GB),
Electron 6 processes RSS 615 MB (PSS 264 MB), total PSS 2.71 GB (2.5–3.0 GB
in six runs). The M0 spike
measured 41.7–46.9 s to a usable window and 1.67 GB PSS without a sysimage:
here the species data, the 19 worker threads and the solver's warm-up load
and compile at the same time as the page, and the compiled solver and the
model stay in memory. With a sysimage (M5) the M0 spike was usable after
2.4 s.

## What remains (M4)

- The tabs: Case (editor), Queue, Results, Compare, Species, with the rest
  of `State` and the plots (on screen StipplePlotly; files from "Save plot"
  with CairoMakie, owner decision PR #13).
- `htg plot` (CairoMakie).
- Roboto as a pinned artifact, like the icon font, so the text looks the same
  on every system (owner decision PR #13).
- A sysimage for a fast start (M5, with packaging), and then the `.theme`
  file of Stipple's precompile step (owner decision PR #13).

## Owner decisions (PR #13)

1. **Plot files: CairoMakie.** On screen the plots are Plotly
   (StipplePlotly); "Save plot" and `htg plot` write PNG, SVG and PDF with
   CairoMakie, sized in mm, from the same plot description (as Python used
   Matplotlib for exports).
2. **`.theme` from Stipple's precompile step:** decided at M5 with the
   sysimage. Until then the install command precompiles from a temporary
   folder and `.theme` is in `.gitignore`.
3. **Warm-up during the page load:** kept; the owner tests start-up again
   with the sysimage.
4. **Immediate exit after "Stopped. ...":** kept.
5. **Text font:** Roboto is bundled later as an artifact, like the icon font.
6. **External requests:** the window keeps blocking all of them; a web link
   a tab needs opens in the system browser.
7. **Status line:** the same text as Python.
