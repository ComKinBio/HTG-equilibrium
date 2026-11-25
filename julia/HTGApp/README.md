# HTGApp.jl

The Julia version of the HTG equilibrium app (`htg gui`; Phase 2, milestone
M4 of `docs/spec/julia_migration_plan.md`). It is built on the science core
`julia/HTGEquilibrium` and on Genie, Stipple and StippleUI (the page, Vue and
Quasar) and Electron.jl (the window), the stack chosen in the M0 spike
(`julia/spikes/m0_window`).

So far it has the shell of the Python app (`src/htg/gui/app.py`): the window,
the header with the title and the status line, the five tabs, the app state
and the app preferences; and the Case tab (the case editor, `gui/editor.py`,
`case_form.py`, `document.py`). The other tabs show a placeholder; their
screens are ported in the next steps of M4. Until then the Python app (`uv run
htg gui`) is the one to use for runs and results.

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
| `src/state.jl` | `State` (model, pool, settings, results and cases folders, `phase_options`, `feed_species`, `document`, `submit`, `comparison_files`), `use_model!`, `load!` (model and pool in the background), `status_text`, `worker_summary`, `output_folder` | `gui/app.py` (`State`), `gui/common.py` |
| `src/case_form.jl` | The Case tab's wording (every caption, tooltip and message), `parse_number`, notices (`conflict_notice`, `entered_sum_notice`, the SUPCRT92 `band_notice`), point counts, the feed check and summary numbers, unit conversions of the feed | `gui/case_form.py` |
| `src/document.jl` | `CaseDocument`: the case being edited (old files converted, dropped solvers, salts, the snapshot for unsaved changes), `open_document`, `new_document`, `save_document` | `gui/document.py` |
| `src/editor.jl` | `CaseEditor`: the whole form without the page (fields, parsing, problems, status line, every operation of the form, `form_view`) | `gui/editor.py` (`CaseEditor`) |
| `src/case_tab.jl` | The Case tab on the page: the Vue/Quasar template, the `case_op` event, one editor per page, native Save as, Check in the background | `gui/editor.py` (the NiceGUI part) |
| `src/page.jl` | The page: header, status line (pushed every second), tabs and panels; the icon font route; the CSS (NiceGUI's colours, the panel width) | `gui/app.py` (`build_page`) |
| `src/window.jl` | Can the window open? (`window_problem`, `explain_electron_failure`, `window_or_browser`) | `gui/window.py` |
| `src/shell.jl` | `run_app`: server, window, life cycle, native file dialogs (for the tabs), the external-request block | `gui/app.py` (`run_app`) |
| `src/cli.jl` | `main` (`julia -m HTGApp gui ...`; the core's `parse_cli`) | `cli.py` (`cmd_gui`) |
| `src/selftest.jl` | The window self-test (`HTG_SELFTEST=1`) | – |
| `test/runtests.jl` | Unit tests (`Pkg.test`) | `tests/test_gui.py` (shell parts), `tests/test_window.py` |
| `test/test_case_tab.jl` | Unit tests of the Case tab (included by `runtests.jl`) | `tests/test_case_tab.py`, `tests/test_case_format.py` (documents) |
| `test/crosscheck_python.jl` | Manual cross-check: every case file saved by the Case tab loads in Python as the same case | – |
| `test/dump_state.py`, `test/data/state_python.json` | Python's `phase_options`, `feed_species` and `worker_summary` texts (7.9 kB, committed) | – |
| `test/selftest.jl` | Runs the app three times with the self-test and checks what is left behind | – |

The app preferences (`settings.toml`) are in the core
(`HTGEquilibrium/src/settings.jl`, see the core README): the same file as the
Python app, so both read each other's choices.

## What is ported

- `State`: `model`, `pool`, `settings`, `results_dir`, `cases_dir` (resolved
  as `Path.resolve` does), `phase_options` and `feed_species` (tested equal to
  Python's), `use_model!`, `submit` (into `results/<case>-<yyyymmdd-HHMMSS>`).
  `document` (the case being edited) and `comparison_files`. Not yet:
  `sources`, `run`, `tables`; they come with the tabs that use them.
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
- **The Case tab**, top to bottom as in Python, with Python's wording: the
  file bar (Case file select with search over the cases folder, New, Save,
  Save as…, the file status, Check, Add to queue), the status line ("Ready to
  run: …" or "✗ <first problem>"), the notices (old format, a file the form
  can't show, a dropped solver, added salts); Name and Notes; Temperature and
  pressure (every T at every P or a list of points, the units as preferences,
  the point count, the collapsible Input formats table whose state is a
  preference); the Feed card in both modes (elemental analysis with kg | %
  for solids and ash, the elements table with Entered and Used, add and
  remove elements, the scaling notice; species amounts with water always
  present, kg | mol with conversion, species selects with search) and the
  live Feed summary; Possible products (chip selects with search per phase);
  Carbon split (free, exactly, at least, at most, between for files that have
  it; solid carbon means); Aqueous phase (with the conflict notice); Solver
  (method, initial guesses with the live tooltip, continuation); Advanced
  settings, collapsed (starting point, solver details, when two results
  differ, minimum amounts, atomic weights, measured data). Every keystroke
  updates the derived numbers; a field that can't be read shows its message
  under it and keeps the case unchanged; Save, Save as and Add to queue are
  disabled while there is a problem, and the status line names the first.
  Check runs the preflight checks of every point in the background with a
  spinner and shows them in a dialog. Add to queue submits through
  `State.submit` and opens the Queue tab.
- `parse_numbers`, `format_numbers` and `parse_points` are in the core
  (`HTGEquilibrium/src/number_lists.jl`), with their Python tests.

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
- **Placeholders.** Each tab but Case says "The <Tab> tab is not in the Julia version
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
- **Case files are written fresh** (owner decision PR #8): saving writes the
  file with the core's `dump_case` (format 2, the core's layout and order);
  comments and the layout of the original file are not kept (Python patched
  the original text with tomlkit, `update_document`, not ported). Added salts
  ([feed.salts]) are read into the case and written back, so they survive a
  save (tested with `cases/manure_ions_minerals.toml`); the notice "This case
  has added salts. The form can't edit them; they are kept when you save."
  stays. A saved file also lists some defaults explicitly (for example
  `model = "electrolyte"`, `method = "ipopt"`); it loads as the same case in
  both apps (`test/crosscheck_python.jl`). Unsaved changes are tracked as in
  Python (a snapshot of the cleaned data).
- **Solvers** (owner decision PR #10): the Method select offers IPOPT
  (default) and SLSQP only; the finite-difference settings and their
  "Finite-difference tolerance/step" fields are gone, and the method tooltip
  is "IPOPT is the default and usually the fastest. Try SLSQP if IPOPT fails
  at some points." A file with `slsqp-fd` or `trust-constr` opens with IPOPT
  and the core's note ("The solver 'slsqp-fd' is no longer available; this
  case now uses 'ipopt' (the default).", and one for `fd_ftol`/`fd_eps`),
  shown with the other notices until the case is saved.
- **Water model** (owner decision PR #10, wording proposed): the Aqueous
  phase card has a "Water model" toggle, SUPCRT92 (the default, no key
  written) | IAPWS-95 (`[aqueous] water_model = "IAPWS-95"`), with the
  tooltip "The water properties of the Electrolyte model: the ion data (HKF)
  and the dielectric constant. SUPCRT92 is the reference and the default.
  IAPWS-95 is the international standard for water and has no failure bands
  near the critical point." It is disabled with the Molecular model ("Used
  only with the Electrolyte model.") and without an aqueous phase ("Used only
  when the case includes an aqueous phase."); the reason is shown under it.
- **SUPCRT92 failure bands** (owner decision PR #10, warn only, wording
  proposed): with the Electrolyte model, SUPCRT92 and an aqueous phase, the
  Temperature and pressure card lists the points in or within 0.5 °C of a
  band at their pressure (`supcrt92_band_check`): "SUPCRT92 can fail at 379,
  386 °C at 240 bar: in or within 0.5 °C of a band where its water density
  does not converge. Results at these points may be wrong. Move these
  temperatures, or choose IAPWS-95 as the water model." (at most 8 points
  listed, then "and N more"; in the chosen units). Points at or above the
  aqueous cut-off are left out (no aqueous phase there). The check runs in a
  background task when the points change (about 0.3 s per new pressure); the
  notice appears when it is done. It does not block saving or queueing.
- **Save as** in the window uses the native save dialog ("Save case file as",
  Case files *.toml; owner decision PR #3). In the browser (no native dialog)
  it shows Python's dialog with the path as text. A name chosen without an
  extension gets `.toml`, so the file appears in the Case file list.
- **Check while the species data load** says "The app is still loading the
  species data. Try again in a moment." (Python's app could not be in that
  state). A species-amounts feed shows the same message under the feed until
  the data are loaded (the molar masses come from them).
- **Unsaved changes survive a page reload**: the document and its snapshot
  live in the app state; a reloaded page (or a second browser tab) gets a new
  form over the same document. Two pages edit the same case; a change in one
  is not pushed to the other until it is reloaded.
- **Notifications** use Quasar's (`$q.notify`) with Python's texts; an error
  message is the exception's message without its type ("Not queued: The app
  is still loading …").
- **How the page talks to the server.** Every change in the form is one
  event to the server, which answers with the new view (labels, messages,
  derived numbers); typing has no extra delay (no 300 ms debounce). The
  events of a page are handled in order in a task of their own, so a slow
  first event (compiling while the solver's warm-up runs) does not hold up
  the page's connection. The form appears a few seconds after the page
  connects; while the solver's warm-up compiles it can take until the
  species data are loaded. It shows "Starting the Case tab." until then.
- **Preferences API.** `settings[key]`, `settings[key] = value`,
  `export_defaults`, `set_export_defaults`; `SettingsKeyError` and
  `SettingsValueError` instead of `KeyError` and `ValueError` (see the core
  README).

## Tests

`julia --project=julia/HTGApp -t 4 -e 'using Pkg; Pkg.test()'`: 290 tests (220 of them the
Case tab's), 1 min 24 s (2026-11-25).

- `worker_summary` against Python's texts; `output_folder`; path resolution.
- `State`: folders, preferences, `phase_options` and `feed_species` equal to
  Python's, the status line with a pool (a stand-in solver), `submit`.
- Loading in the background: good data, and a broken data folder.
- `gui` options, usage and help; the window check and its messages.
- The page as served by Genie: no reference to the internet, no Google
  font, no Stipple theme files, the font route, and no `.theme` file written
  into the working folder.
- The Case tab (`test_case_tab.jl`, the tests of Python's `test_case_tab.py`
  and more): the WP4 feed numbers (10 kg, 15 %, 8.45 % ash, C 30.4 H 4.41 O
  53.58 N 2.66 S 0.51: Used 33.20, 4.82, 58.52, 2.91, 0.56; water 8.500 kg /
  471.8 mol; ash 0.127 kg; ash-free solids 1.373 kg; C 0.456 kg / 37.96 mol;
  the 91.56 % notice); the feed errors and their fields; the conflict notice
  (its three conditions, in K); the initial-guesses tooltip (1 guess, 5
  guesses); point counts, element symbols, `parse_number`, `file_status`; the
  SUPCRT92 band notice (which cases are checked; the real check at 240 bar:
  379 °C in a band, 386 °C near one; the wording, units, the list limit);
  documents (every example case opens and saves back to an equal case, salts
  kept, old files converted, dropped solvers, a file the form can't show);
  the editor driven as the page drives it: the feed card, ash in kg and the
  ash error blocking Save and Add to queue, elements added and removed, the
  conflict notice following each condition and the temperature unit, the
  initial guesses, conditions and point counts in both modes, the units and
  Input formats as preferences, the band notice (in the background too), the
  water model and when it can be chosen, the simple selects, species amounts
  (kg | mol, repeated or missing species, switching modes and back), open,
  Save, Save as, a reloaded page, salts and blocked files, Check, Add to
  queue (while loading and with a stand-in pool); the page's markup.

**Cross-check with Python** (manual; needs uv):
`julia --project=julia/HTGApp/test julia/HTGApp/test/crosscheck_python.jl`
opens every case file of the app (cases/, tests/data/legacy_cases, the
template) in a CaseEditor, changes a value and changes it back, saves it with
Save as, and has Python load the saved file and the original and compare the
two cases as Python's `dump_case` writes them (the name and the dropped
solver's `fd_ftol`/`fd_eps` left out). Result (2026-11-25): all 12 the same.

**Self-test of the window** (needs a desktop session; Linux for the process
checks):

```
julia julia/HTGApp/test/selftest.jl              # three sessions, ~8 min
julia julia/HTGApp/test/selftest.jl --offline    # the same without network (unshare -rn)
```

It starts the app three times as `htg gui` does, with `HTG_SELFTEST=1`, one
new config folder and an empty working folder each time (and a temporary
results folder). In the window it checks: the title, the
fixed header (title left, status line right), the header colour, the tab
names and icons, that the icons are glyphs of the bundled font, that the tabs
switch (the Case tab shows the form, the others the placeholder), the panel width, that every resource came from the
app's server and that no request was blocked; that the data load and the pool
start in the background, and the status line while loading; that the status
line follows the pool (`set_workers`); that the warm-up finishes. The Case
tab (sessions 1 and 2): the WP4 feed typed into the fields (the Used column,
the entered sum and its notice, the Feed summary); ash switched to kg (the
value converted) and set to 2 kg (the message under the field, Add to queue,
Save and Save as disabled, the status line); an example case opened,
Initial guesses changed (unsaved changes), Save as (the native dialog is
answered by the test: its title and suggested file are checked) and the file
checked; the page reloaded (the saved case still in the form); a new case
added to the queue (the Queue tab opens, the run finishes). The first
session sets the pressure unit to MPa; the second must read it at start-up.
The third closes the window during the warm-up; its process must end within
seconds. After each session: exit code 0, no Julia or Electron process left,
nothing written to the working folder, no error or warning logged, and no
`htg-*` folder left in the system temp folder.

Result (2026-11-25, Linux, Wayland): 34/34, 33/33 and 11/11 checks passed
(the setting was kept; the queued case ran in both sessions), no process,
file or log message left, and the third process ended 0.2 s after its window
closed. The `--offline` run (no network at all, loopback only) was last done
on 2026-11-24, before the Case tab: the same checks passed then, the page had
14 resources, all from the app's server, and nothing was blocked.

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
| Julia running, HTGApp loaded | 6.4 s | 6.8 s |
| Window open with the start-up page | 14.1 s | 15.2 s |
| Window usable (page loaded, websocket connected) | 56.1 s | 63.0 s |
| Species data loaded, pool started (as seen by the self-test) | 75.2 s | 63.8 s |
| Warm-up done (solver compiled, "Workers: 19 idle") | 174 s | 188 s |

(2026-11-25, with the Case tab; before it, on 2026-11-24: usable after 62.3
and 56.1 s, warm-up done after 179 and 176 s.) Idle memory after the warm-up
(session 2, after the Case tab's checks and one queued run): Julia RSS
2.88 GB (PSS 2.79 GB), Electron 6 processes RSS 676 MB (PSS 316 MB), total
PSS 3.10 GB (2.71 GB on 2026-11-24 without the Case tab). The M0 spike
measured 41.7–46.9 s to a usable window and 1.67 GB PSS without a sysimage:
here the species data, the 19 worker threads and the solver's warm-up load
and compile at the same time as the page, and the compiled solver and the
model stay in memory. With a sysimage (M5) the M0 spike was usable after
2.4 s.

## What remains (M4)

- The Case tab: the Species tab's "Add to case" (Python's `add_feed_species`
  and adding products to the open case) comes with the Species tab; the
  `on_queued` hook of Python's editor (the Queue tab selecting the new run)
  with the Queue tab.
- The tabs: Queue, Results, Compare, Species, with the rest
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
