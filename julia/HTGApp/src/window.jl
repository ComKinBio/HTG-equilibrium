# Can the app open in its own window (Electron)? (Python: gui/window.py, which
# checked pywebview and Qt.)
#
# Checked before the window opens, so that the app can fall back to the browser
# with one line saying why. On Linux a missing display or a missing system
# library makes Electron stop at once; a tiny Electron program is started first
# (about 0.2 s) and its error output explained. On Windows and macOS only the
# Electron binary is checked (Python made no check there either); an Electron
# that still fails to start is caught when the window opens.

const ELECTRON_LINUX_LIBRARIES = "libnss3, libgtk-3-0, libgbm1 and libasound2"

const PROBE_JS = """
const { app } = require('electron');
app.whenReady().then(() => app.exit(0));
"""

"""The Electron binary Electron.jl uses, or nothing when it is not installed."""
function electron_binary()
    path = try
        Electron.get_electron_binary_cmd()
    catch
        return nothing
    end
    return isfile(path) ? path : nothing
end

"""Nothing if the window can open; otherwise one line saying what is missing and how to fix it."""
function window_problem(; timeout=30)
    binary = electron_binary()
    binary === nothing &&
        return "Electron is not installed. Run: julia --project=julia/HTGApp -e \"using Pkg; Pkg.instantiate()\"."
    Sys.islinux() || return nothing
    if isempty(get(ENV, "DISPLAY", "")) && isempty(get(ENV, "WAYLAND_DISPLAY", ""))
        return "no display was found. Run htg gui from a desktop session."
    end
    dir = mktempdir(; prefix="htg-probe-")
    try
        script = joinpath(dir, "probe.js")
        write(script, PROBE_JS)
        err = Pipe()
        cmd = pipeline(ignorestatus(`$binary --no-sandbox $script`); stdout=devnull, stderr=err)
        proc = run(cmd; wait=false)
        close(err.in)
        text = @async read(err, String)
        if timedwait(() -> process_exited(proc), Float64(timeout); pollint=0.02) !== :ok
            kill(proc)
            return "Electron did not start within $timeout s."
        end
        success(proc) && return nothing
        return explain_electron_failure(fetch(text))
    catch exc
        exc isa InterruptException && rethrow()
        return "Electron could not be checked ($(sprint(showerror, exc)))."
    finally
        rm(dir; recursive=true, force=true)
    end
end

"""One line for the error output of a failed Electron start (Linux)."""
function explain_electron_failure(stderr_text)
    text = something(stderr_text, "")
    if occursin("Missing X server", text) || occursin("cannot open display", text) ||
       occursin("Failed to connect to Wayland", text) || occursin("platform failed to initialize", text)
        return "no display was found. Run htg gui from a desktop session."
    end
    library = match(r"(lib[\w+.-]+?\.so[\d.]*): cannot open shared object file", text)
    if library !== nothing
        return "the system library $(library[1]) is missing. Install it with your package manager " *
               "(Electron needs, among others, $ELECTRON_LINUX_LIBRARIES on Debian or Ubuntu)."
    end
    lines = [strip(l) for l in split(text, '\n') if !isempty(strip(l))]
    reason = isempty(lines) ? "Electron stopped without a message" : first(lines[end], 150)
    return "Electron could not start ($reason)."
end

"""
    window_or_browser(browser=false; out=stdout) -> Bool

True to open the app's own window. Prints one line and returns false when it
cannot open (Python: `native_or_browser`).
"""
function window_or_browser(browser=false; out::IO=stdout, problem=() -> window_problem())
    browser && return false
    reason = problem()
    reason === nothing && return true
    println(out, "Cannot open a separate window: $reason Opening the app in the browser instead " *
                 "(htg gui --browser does this directly).")
    return false
end
