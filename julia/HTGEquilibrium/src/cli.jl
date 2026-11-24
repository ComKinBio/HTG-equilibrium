# Command line interface (Python: htg/cli.py).
#
#     htg run CASE.toml [CASE2.toml ...] [--out DIR] [--workers N] [--no-plots]
#     htg check CASE.toml
#     htg species [--all] [--phase PHASE] [--search TEXT]
#     htg new CASE.toml
#     htg plot RESULTS_DIR          (not in the Julia version yet: exits with a message)
#     htg gui [...]                 (the app; in the package julia/HTGApp, see below)
#
# run, check and species read your own species database too (see
# user_species.jl); --user-data DIR chooses its folder, --no-user-species turns
# it off. `--data DIR` (before the command) chooses the database folder.
#
# Started with `julia/htg` (Linux, macOS) or `julia\htg.cmd` (Windows), which run
# `julia --project=<this package> --threads=auto -m HTGEquilibrium ARGS...`;
# for the `gui` command they run `-m HTGApp` with the project julia/HTGApp
# instead, so that this package does not load the GUI packages (HTGApp parses
# the arguments with `parse_cli` below). Started here, `gui` exits with a message
# that says how to start the app.
# `main(args)` is the entry point (`@main`). Arguments, defaults, printed text and
# exit codes are Python's (argparse: usage errors exit 2 with "htg: error: ...";
# ValueError and missing files print "Error: ..." and exit 2; failed points or
# preflight problems exit 1). `--workers N` means N worker threads; Julia's thread
# count is fixed at start-up, so N is capped at the threads Julia has (the
# launchers start Julia with one thread per CPU thread).

"""An exit of the command line with a code (argparse's `parser.exit`)."""
struct _CLIExit <: Exception
    code::Int
end

"""A Python FileNotFoundError ("[Errno 2] No such file or directory: 'x'")."""
struct _FileNotFound <: Exception
    path::String
end
Base.showerror(io::IO, e::_FileNotFound) = print(io, "[Errno 2] No such file or directory: ", py_repr(e.path))

# ---------------------------------------------------------------- argument parsing (argparse subset)
struct _Opt
    flags::Vector{String}    # e.g. ["-h", "--help"]
    dest::String
    kind::Symbol             # :help, :flag (store_true), :value, :int
    default::Any             # flags: false
end
_Opt(flags, dest, kind) = _Opt(flags, dest, kind, kind == :flag ? false : nothing)

struct _Command
    name::String
    usage::String            # argparse's usage text (also printed with errors)
    help::String             # the text after the usage for -h
    options::Vector{_Opt}
    positionals::Vector{Pair{String,Symbol}}   # name => :one or :plus
end

const _HELP_OPT = _Opt(["-h", "--help"], "help", :help)
const _USER_OPTS = [_Opt(["--user-data"], "user_data", :value), _Opt(["--no-user-species"], "no_user_species", :flag)]
const _USER_HELP = """
  --user-data DIR    folder of your own species database (default:
                     \$HTG_USER_DATA, else the htg-equilibrium user data
                     folder)
  --no-user-species  use the shipped databases only, without your own species
"""

const _MAIN_USAGE = "usage: htg [-h] [--data DATA] {run,check,species,new,plot,gui} ..."
const _MAIN = _Command("htg", _MAIN_USAGE, """

Hydrothermal gasification equilibrium.

positional arguments:
  {run,check,species,new,plot,gui}
    run                 solve case files and write tables (and plots: not in
                        the Julia version yet)
    check               fast checks of a case file (no equation of state)
    species             list species that can be used in case files
    new                 write a commented case-file template
    plot                re-plot saved results (not in the Julia version yet)
    gui                 open the app in its own window (its tabs are not in
                        the Julia version yet)

options:
  -h, --help            show this help message and exit
  --data DATA           database folder (default: the repository data/ folder)
""", [_HELP_OPT, _Opt(["--data"], "data", :value)], ["command" => :one])

const _COMMANDS = OrderedDict{String,_Command}(
    "run" => _Command("htg run", """
usage: htg run [-h] [--out OUT] [--workers WORKERS] [--no-plots]
               [--user-data DIR] [--no-user-species]
               case [case ...]""", """

positional arguments:
  case

options:
  -h, --help         show this help message and exit
  --out OUT          output folder (one case), or parent folder (several
                     cases); default: results/<case name>
  --workers WORKERS  worker threads; independent points and continuation
                     chains of all cases are shared between them (at most
                     the threads Julia was started with)
  --no-plots         do not write plots (the Julia version writes none yet)
""" * _USER_HELP,
        [_HELP_OPT, _Opt(["--out"], "out", :value), _Opt(["--workers"], "workers", :int, 1),
         _Opt(["--no-plots"], "no_plots", :flag), _USER_OPTS...], ["case" => :plus]),
    "check" => _Command("htg check", "usage: htg check [-h] [--user-data DIR] [--no-user-species] case", """

positional arguments:
  case

options:
  -h, --help         show this help message and exit
""" * _USER_HELP, [_HELP_OPT, _USER_OPTS...], ["case" => :one]),
    "species" => _Command("htg species", """
usage: htg species [-h] [--all] [--phase PHASE] [--search SEARCH]
                   [--user-data DIR] [--no-user-species]""", """

options:
  -h, --help         show this help message and exit
  --all              include species with incomplete data
  --phase PHASE      e.g. molecular, organic, solid, ion
  --search SEARCH    part of a name or formula
""" * _USER_HELP,
        [_HELP_OPT, _Opt(["--all"], "all", :flag), _Opt(["--phase"], "phase", :value),
         _Opt(["--search"], "search", :value), _USER_OPTS...], Pair{String,Symbol}[]),
    "new" => _Command("htg new", "usage: htg new [-h] case", """

positional arguments:
  case

options:
  -h, --help  show this help message and exit
""", [_HELP_OPT], ["case" => :one]),
    "plot" => _Command("htg plot", "usage: htg plot [-h] results", """

Not available in the Julia version yet (use the Python app: uv run htg plot).

positional arguments:
  results

options:
  -h, --help  show this help message and exit
""", [_HELP_OPT], ["results" => :one]),
    "gui" => _Command("htg gui", """
usage: htg gui [-h] [--browser] [--workers WORKERS] [--port PORT]
               [--host HOST] [--results RESULTS] [--cases CASES]
               [--no-browser] [--user-data DIR] [--no-user-species]""", """

The Julia version has the window, header, status line and tabs; the tabs
themselves are not ported yet (for them use the Python app: uv run htg gui).

options:
  -h, --help         show this help message and exit
  --browser          open the app in the web browser instead
  --workers WORKERS  worker threads (default: from free cores and memory; at
                     most the threads Julia was started with)
  --port PORT        default: 8080 in the browser, a free port for the window
  --host HOST        network address (default 127.0.0.1: only this computer;
                     0.0.0.0: everyone on the network)
  --results RESULTS  folder for results (default: results)
  --cases CASES      folder with case files (default: cases)
  --no-browser       only start the server; open the address in a browser
                     yourself
""" * _USER_HELP,
        [_HELP_OPT, _Opt(["--browser"], "browser", :flag), _Opt(["--workers"], "workers", :int),
         _Opt(["--port"], "port", :int), _Opt(["--host"], "host", :value, "127.0.0.1"),
         _Opt(["--results"], "results", :value, "results"), _Opt(["--cases"], "cases", :value, "cases"),
         _Opt(["--no-browser"], "no_browser", :flag), _USER_OPTS...],
        Pair{String,Symbol}[]),
)

"""argparse's error: usage and "prog: error: message" on stderr, exit 2."""
function _usage_error(cmd::_Command, message, err::IO)
    println(err, cmd.usage)
    println(err, "$(cmd.name): error: $message")
    throw(_CLIExit(2))
end

const _NEGATIVE_NUMBER = r"^-\d+$|^-\d*\.\d+$"

"""
How argparse reads `arg` with `cmd`'s options: (:pos,), (:opt, option, explicit
value or nothing) or (:unknown,). Long options may be abbreviated (unique prefix).
"""
function _classify(cmd::_Command, arg::AbstractString, err::IO)
    (isempty(arg) || arg[1] != '-') && return (:pos, nothing, nothing)
    for o in cmd.options
        arg in o.flags && return (:opt, o, nothing)
    end
    length(arg) == 1 && return (:pos, nothing, nothing)
    name, explicit = occursin('=', arg) ? split(arg, '='; limit=2) : (arg, nothing)
    if explicit !== nothing
        for o in cmd.options
            name in o.flags && return (:opt, o, String(explicit))
        end
    end
    if startswith(name, "--")
        matches = [(o, f) for o in cmd.options for f in o.flags if startswith(f, "--") && startswith(f, name)]
        if length(matches) > 1
            _usage_error(cmd, "ambiguous option: $name could match $(join(last.(matches), ", "))", err)
        elseif length(matches) == 1
            return (:opt, matches[1][1], explicit === nothing ? nothing : String(explicit))
        end
    end
    occursin(_NEGATIVE_NUMBER, arg) && return (:pos, nothing, nothing)
    occursin(' ', arg) && return (:pos, nothing, nothing)
    return (:unknown, nothing, nothing)
end

"""Python's `int(text)` for an option value; nothing if it is not an integer."""
function _parse_int(text)
    s = replace(strip(text), r"(?<=\d)_(?=\d)" => "")
    return tryparse(Int, s)
end

"""
Parse `args` with `cmd` (argparse rules for these commands: options anywhere,
`--opt=value`, abbreviations, positionals in runs between options, `--`).
Returns (values, unrecognized arguments). -h prints the help and exits 0.
"""
function _parse(cmd::_Command, args, out::IO, err::IO)
    values = Dict{String,Any}(o.dest => o.default for o in cmd.options if o.kind != :help)
    # After the first "--" every word is a positional. The "--" itself (:sep) is
    # dropped from a positional's values but kept among unrecognized arguments
    # (argparse 3.12).
    words, kinds = String[], Any[]
    plain = false
    for a in args
        push!(words, a)
        if !plain && a == "--"
            plain = true
            push!(kinds, (:sep, nothing, nothing))
        else
            push!(kinds, plain ? (:pos, nothing, nothing) : _classify(cmd, a, err))
        end
    end
    extras = String[]
    next_positional = 1
    i = 1
    while i <= length(words)
        kind, opt, explicit = kinds[i]
        if kind in (:pos, :sep)
            j = i
            while j <= length(words) && kinds[j][1] in (:pos, :sep)
                j += 1
            end
            group, seps = words[i:j-1], [k[1] == :sep for k in kinds[i:j-1]]
            k = 1
            while next_positional <= length(cmd.positionals) && any(!s for s in seps[k:end])
                name, nargs = cmd.positionals[next_positional]
                if nargs == :plus
                    values[name] = [w for (w, s) in zip(group[k:end], seps[k:end]) if !s]
                    k = length(group) + 1
                else
                    while seps[k]   # "--" before the value
                        k += 1
                    end
                    values[name] = group[k]
                    k += 1
                    while k <= length(group) && seps[k]
                        k += 1
                    end
                end
                next_positional += 1
            end
            append!(extras, group[k:end])
            i = j
        elseif kind == :unknown
            push!(extras, words[i])
            i += 1
        else
            i = _take_option!(values, cmd, opt, explicit, words, kinds, i, out, err)
        end
    end
    missing_names = [name for (name, _) in cmd.positionals[next_positional:end]]
    isempty(missing_names) ||
        _usage_error(cmd, "the following arguments are required: $(join(missing_names, ", "))", err)
    return values, extras
end

"""Store option `opt` found at word `i`; returns the index of the next word."""
function _take_option!(values, cmd, opt, explicit, words, kinds, i, out, err)
    label = join(opt.flags, "/")
    if opt.kind == :help
        print(out, cmd.usage, "\n", cmd.help)
        throw(_CLIExit(0))
    elseif opt.kind == :flag
        explicit === nothing ||
            _usage_error(cmd, "argument $label: ignored explicit argument $(py_repr(explicit))", err)
        values[opt.dest] = true
        return i + 1
    end
    if explicit !== nothing
        value = explicit
    elseif i < length(words) && kinds[i+1][1] == :pos
        value = words[i+1]
        i += 1
    else
        _usage_error(cmd, "argument $label: expected one argument", err)
    end
    if opt.kind == :int
        number = _parse_int(value)
        number === nothing && _usage_error(cmd, "argument $label: invalid int value: $(py_repr(value))", err)
        values[opt.dest] = number
    else
        values[opt.dest] = value
    end
    return i + 1
end

"""
The main parser: options before the command, the command, and the words after
it (which belong to the command's parser, as argparse's subparsers).
"""
function _parse_main(args, out::IO, err::IO)
    values = Dict{String,Any}("data" => nothing)
    separator = something(findfirst(==("--"), args), length(args) + 1)
    kinds = Any[i < separator ? _classify(_MAIN, a, err) : (:pos, nothing, nothing) for (i, a) in enumerate(args)]
    extras = String[]
    i = 1
    while i <= length(args)
        kind, opt, explicit = kinds[i]
        if args[i] == "--"   # argparse: the next word is the command
            i += 1
            i <= length(args) && return args[i], args[i+1:end], values, extras
        elseif kind == :pos
            return args[i], args[i+1:end], values, extras
        elseif kind == :unknown
            push!(extras, args[i])
            i += 1
        else
            i = _take_option!(values, _MAIN, opt, explicit, args, kinds, i, out, err)
        end
    end
    _usage_error(_MAIN, "the following arguments are required: command", err)
end

"""Parse the whole command line: (command name, values)."""
function parse_cli(args, out::IO=stdout, err::IO=stderr)
    command, rest, main_values, extras = _parse_main(collect(String, args), out, err)
    haskey(_COMMANDS, command) ||
        _usage_error(_MAIN, "argument command: invalid choice: $(py_repr(command)) " *
                            "(choose from $(join(keys(_COMMANDS), ", ")))", err)
    values, sub_extras = _parse(_COMMANDS[command], rest, out, err)
    append!(extras, sub_extras)
    isempty(extras) || _usage_error(_MAIN, "unrecognized arguments: $(join(extras, " "))", err)
    values["data"] = main_values["data"]
    return command, values
end

# ---------------------------------------------------------------- helpers
"""`str(pathlib.Path(text))`: repeated and trailing separators and `.` parts dropped."""
function _pypath(text::AbstractString)
    sep = Sys.iswindows() ? "\\" : "/"
    text = Sys.iswindows() ? replace(text, "/" => "\\") : text
    isempty(text) && return "."
    absolute = startswith(text, sep)
    parts = [p for p in split(text, sep) if !isempty(p) && p != "."]
    path = join(parts, sep)
    absolute && (path = sep * path)
    return isempty(path) ? "." : path
end

_pyjoin(a, b) = _pypath(a) == "." ? _pypath(b) : _pypath(a) * (Sys.iswindows() ? "\\" : "/") * _pypath(b)

"""The user species folder for this command (nothing with --no-user-species)."""
_user_data(values) = resolve_user_dir(get(values, "user_data", nothing); enabled=!get(values, "no_user_species", false))

_model(values) = load_model(values["data"]; user_data=_user_data(values))

function _load_case(path)
    isfile(path) || throw(_FileNotFound(_pypath(path)))
    return load_case(path)
end

"""Python's `_print_table`: pandas `to_string(index=False)` of the table (or of `columns`)."""
_print_table(out::IO, table::DataTable, columns=nothing) =
    println(out, to_string(columns === nothing ? table : select(table, columns)))

"""The preflight table as pandas has it: `aqueous_removed` is None where it was not set."""
function _check_table(table::DataTable)
    table = copy(table)
    if haskey(table, "aqueous_removed")
        table.data["aqueous_removed"] = Any[v === missing ? nothing : v for v in table["aqueous_removed"]]
    end
    return table
end

"""The number of worker threads for `--workers N`: N, at most the threads Julia has."""
function worker_threads(requested::Integer, out::IO=stdout)
    available = Threads.nthreads(:default)
    requested <= available && return Int(requested)
    println(out, "Julia has $available thread(s); using $available worker(s) instead of $requested " *
                 "(start Julia with more threads: julia -t $requested, or the htg launcher, which uses all CPU threads)")
    return available
end

# ---------------------------------------------------------------- commands
function cmd_run(values, out::IO, err::IO)
    cases = [_load_case(path) for path in values["case"]]
    names = [c.name for c in cases]
    length(Set(names)) == length(names) ||
        throw(CaseError("Case names must differ (they name the output folders): $(py_repr(names))"))
    out_arg = values["out"]
    outputs = if out_arg !== nothing && !isempty(out_arg) && length(cases) == 1
        [_pypath(out_arg)]
    else
        [_pyjoin(out_arg === nothing || isempty(out_arg) ? "results" : out_arg, c.name) for c in cases]
    end
    values["no_plots"] ||
        println(out, "Plots are not written by the Julia version yet (as with --no-plots); " *
                     "for plots use the Python app: uv run htg plot RESULTS_DIR")
    model = _model(values)
    for case in cases
        println(out, "Case $(case.name): $(length(grid(case))) point(s), $(length(chains(case))) independent piece(s), " *
                     "solver $(case.solver)")
        table = _check_table(check(model, case))
        invalid = take_rows(table, [s != "ready" for s in table["status"]])
        if nrows(invalid) > 0
            println(out, "Preflight found problems at some points (they will probably fail):")
            _print_table(out, invalid, ["T_C", "P_bar", "error"])
        end
    end
    start = time()
    workers = values["workers"] > 1 ? worker_threads(values["workers"], out) : 1
    if workers > 1
        function on_point(run, T, P, outcome, start_label)
            println(out, "$(run.name): ", point_text(T, P, outcome, start_label))
            flush(out)
        end
        println(out, "Starting $workers worker threads (the first results take about a minute while Julia compiles)")
        sweeps = WorkerPool(workers; model, warm_up=false, on_point) do pool
            runs = [submit(pool, case, output) for (case, output) in zip(cases, outputs)]
            foreach(wait, runs)
            [sweep_result(run) for run in runs]
        end
    else
        sweeps = SweepResult[]
        for (case, output) in zip(cases, outputs)
            result = sweep(model, case; progress=line -> (println(out, line); flush(out)))
            save(result, output)
            push!(sweeps, result)
        end
    end
    println(out, "Finished in $(@sprintf("%.0f", time() - start)) s")
    ok = true
    for (case, result, output) in zip(cases, sweeps, outputs)
        table = summary(result)
        if haskey(table, "error")
            setcolumn!(table, "error", [_isnull(v) ? "" : first(String(v), 90) for v in table["error"]])
        end
        println(out, "\n$(case.name):")
        columns = [c for c in ("T_C", "P_bar", "status", "gibbs_J_per_mol_C", "gas_carbon_percent", "aqueous_removed",
                               "error") if haskey(table, c)]
        _print_table(out, table, columns)
        println(out, "Results written to $output/")
        ok = ok && all(==("solved"), table["status"])
    end
    return ok ? 0 : 1
end

function cmd_check(values, out::IO, err::IO)
    case = _load_case(values["case"])
    table = _check_table(check(_model(values), case))
    _print_table(out, table)
    ok = all(==("ready"), table["status"])
    println(out, ok ? "All temperatures passed the preflight checks." :
                      "Some temperatures are invalid; see the error column.")
    return ok ? 0 : 1
end

const _CATALOGUE_COLUMNS = ["Species", "Formula", "Phase", "Ready", "Reason"]

function cmd_species(values, out::IO, err::IO)
    phase = values["phase"]
    if phase !== nothing && !isempty(phase)
        try
            Regex(phase, "i")
        catch exc
            exc isa ErrorException || rethrow()
            throw(CaseError("--phase $(py_repr(phase)) is not a valid regular expression: $(exc.msg)"))
        end
    end
    rows = species(_model(values); include_unavailable=values["all"], phase, search=values["search"])
    table = DataTable([OrderedDict{String,Any}(c => getfield(r, Symbol(c)) for c in _CATALOGUE_COLUMNS) for r in rows];
                      columns=_CATALOGUE_COLUMNS)
    _print_table(out, table, values["all"] ? _CATALOGUE_COLUMNS : _CATALOGUE_COLUMNS[1:3])
    println(out, "$(length(rows)) species")
    return 0
end

function cmd_new(values, out::IO, err::IO)
    path = _pypath(values["case"])
    if ispath(path)
        println(err, "$path already exists; choose another name.")
        return 1
    end
    folder = dirname(path)
    isempty(folder) || isdir(folder) || throw(_FileNotFound(path))
    stem = splitext(basename(path))[1]
    text = replace(read(CASE_TEMPLATE, String), "name = \"my_case\"" => "name = \"$stem\"")
    write(path, text)
    println(out, "Wrote $path. Edit it, then: htg check $path  and  htg run $path")
    return 0
end

function _not_available(command, example)
    function (values, out::IO, err::IO)
        println(err, "Error: htg $command is not available in the Julia version yet (it comes with the plots " *
                     "and the app). For now use the Python version: $example")
        return 2
    end
end

"""`htg gui` started with this package: the app is the package julia/HTGApp (it loads the GUI packages)."""
function cmd_gui(values, out::IO, err::IO)
    println(err, "Error: the app is in the package julia/HTGApp. Start it with the launcher: julia/htg gui " *
                 "(Windows: julia\\htg.cmd gui), or: julia --project=julia/HTGApp -t auto -m HTGApp gui")
    return 2
end

const _HANDLERS = Dict{String,Function}(
    "run" => cmd_run, "check" => cmd_check, "species" => cmd_species, "new" => cmd_new,
    "plot" => _not_available("plot", "uv run htg plot RESULTS_DIR"), "gui" => cmd_gui)

"""Python's `str(exc)` of the errors the command line reports as "Error: ..."."""
_error_text(e::ArgumentError) = e.msg
_error_text(e) = sprint(showerror, e)

"""Errors printed as "Error: ..." with exit code 2 (Python: ValueError, FileNotFoundError)."""
_reported(e) = e isa ValueError || e isa ArgumentError || e isa _FileNotFound || e isa TOML.ParserError

"""
    main(args=ARGS; out=stdout, err=stderr) -> exit code

The `htg` command line (Python: `htg.cli.main`). Run as
`julia --project=julia/HTGEquilibrium -t auto -m HTGEquilibrium run case.toml`
or with the launchers `julia/htg` and `julia/htg.cmd`.
"""
function main(args=ARGS; out::IO=stdout, err::IO=stderr)
    try
        command, values = parse_cli(collect(String, args), out, err)
        try
            return _HANDLERS[command](values, out, err)
        catch exc
            _reported(exc) || rethrow()
            println(err, "Error: ", _error_text(exc))
            return 2
        end
    catch exc
        exc isa _CLIExit && return exc.code
        exc isa InterruptException && rethrow()
        # An unexpected error (Python: a traceback, exit code 1)
        showerror(err, exc, catch_backtrace())
        println(err)
        return 1
    end
end

@main
