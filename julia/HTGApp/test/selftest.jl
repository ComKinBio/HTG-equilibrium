# Self-test of the app window (needs a desktop session; Linux for the process
# and memory checks). Not part of Pkg.test().
#
#   julia julia/HTGApp/test/selftest.jl [--offline] [--idle S]
#
# Starts the app twice as `htg gui` does (`julia -t auto -m HTGApp gui
# --no-user-species`) with HTG_SELFTEST=1, one new config folder for both and
# an empty working folder each. In the window it checks the page
# (src/selftest.jl), the Case tab included (its queued case writes its
# results to a temporary folder given with --results); the first session writes a setting (pressure unit MPa),
# the second must read it at start-up. After each session: exit code 0, every
# check passed, no Julia or Electron process of the session left, nothing
# written to the working folder, no error or warning logged, no htg-* folder
# (session files, window check) left in the system temp folder. A third
# session closes the window during the warm-up; the process must end within
# seconds.
#
# --offline  run the sessions in a network namespace without internet
#            (`unshare -rn`, loopback only), so the page cannot reach the network
# --idle S   in the second session, sit idle S seconds after the warm-up, then
#            report the memory of all processes (default 10)
#
# Prints the time from launch to a usable window (page loaded, websocket
# connected, a server handler ran) and the idle memory.

const APP = normpath(joinpath(@__DIR__, ".."))

function session(n, config; set="", expect="", idle=0, offline=false, early=false)
    work = mktempdir()
    results = mktempdir()   # the Case tab queues a case: its results go here, not to the working folder
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    cmd = `$julia --project=$APP --threads=auto --startup-file=no -m HTGApp gui --no-user-species --results $results`
    offline && (cmd = `unshare -rn sh -c 'ip link set lo up && exec "$@"' sh $cmd`)
    env = copy(ENV)
    merge!(env, Dict("HTG_SELFTEST" => "1", "HTG_SELFTEST_SET" => set, "HTG_SELFTEST_EXPECT" => expect,
                     "HTG_SELFTEST_IDLE" => string(idle), "HTG_SELFTEST_EARLY" => early ? "1" : "", "HTG_CONFIG_DIR" => config, "HTG_T0" => string(time())))
    temp_before = Set(readdir(tempdir()))
    output = IOBuffer()
    t0 = time()
    proc = run(pipeline(ignorestatus(setenv(Cmd(cmd; dir=work), env)); stdout=output, stderr=output))
    elapsed = time() - t0
    text = String(take!(output))
    sleep(1)
    lines = split(text, '\n')
    for line in lines
        occursin(r"^(\[htg\]|Stopped|Starting|Error|Cannot|┌|│|└)", line) && println("  ", line)
    end
    pids = Int[]
    for line in lines
        m = match(r"^\[htg\] processes: julia (\d+); electron ([\d ]+)$", line)
        m === nothing || append!(pids, parse.(Int, vcat(m[1], split(m[2]))))
    end
    left = [p for p in pids if isdir("/proc/$p")]
    usable = match(r"([\d.]+) s  USABLE", text)
    memory = match(r"memory: (.*)", text)
    checks = count(l -> occursin(r"^\[htg\]   (PASS|FAIL)", l), lines)
    failed = count(l -> occursin(r"^\[htg\]   FAIL", l), lines)
    written = readdir(work)
    logged = count(l -> occursin(r"^┌ (Error|Warning)", l), lines)
    temp_left = [d for d in setdiff(readdir(tempdir()), temp_before) if startswith(d, "htg-")]
    # closed during the warm-up: the process ends within seconds of the window (not when the compilation is done)
    closed = match(r"([\d.]+) s  electron exited", text)
    after_close = closed === nothing ? Inf : elapsed - parse(Float64, closed[1])
    quick = !early || after_close < 5
    ok = proc.exitcode == 0 && checks > 0 && failed == 0 && !isempty(pids) && isempty(left) && isempty(written) &&
         logged == 0 && isempty(temp_left) && quick
    println("session $n: exit code $(proc.exitcode); $(checks - failed)/$checks checks passed; " *
            "processes left: $(isempty(pids) ? "not reported" : length(left)) of $(length(pids)); " *
            "files in the working folder: $(isempty(written) ? "none" : join(written, ", ")); " *
            "logged errors and warnings: $logged; " *
            "htg folders left in $(tempdir()): $(isempty(temp_left) ? "none" : join(temp_left, ", ")); " *
            "usable after $(usable === nothing ? "?" : usable[1]) s; process ended after $(round(elapsed; digits=1)) s, " *
            "$(round(after_close; digits=1)) s after the window" * (memory === nothing ? "" : "; " * memory[1]))
    rm(work; recursive=true, force=true)
    rm(results; recursive=true, force=true)
    return ok
end

function main(args)
    offline = "--offline" in args
    i = findfirst(==("--idle"), args)
    idle = i === nothing ? 10 : parse(Float64, args[i+1])
    println("Self-test of the app window", offline ? " (no network)" : "")
    config = mktempdir()
    ok1 = session(1, config; set="MPa", expect="bar", offline)
    ok2 = session(2, config; expect="MPa", idle, offline)
    ok3 = session(3, config; expect="MPa", offline, early=true)
    stored = readdir(config) == ["settings.toml"]
    stored || println("config folder should hold settings.toml only: $(readdir(config))")
    passed = ok1 && ok2 && ok3 && stored
    println(passed ? "Self-test passed." : "Self-test FAILED.")
    return passed ? 0 : 1
end

exit(main(ARGS))
