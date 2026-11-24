# Memory of the current process from /proc (Linux): RSS, PSS (shared pages split
# between the processes that map them) and peak RSS, in MB.
function memory_report(pid = getpid())
    status = read("/proc/$pid/status", String)
    kb(key, text) = (m = match(Regex("^$key:\\s+(\\d+) kB", "m"), text); m === nothing ? 0 : parse(Int, m[1]))
    rollup = isfile("/proc/$pid/smaps_rollup") ? read("/proc/$pid/smaps_rollup", String) : ""
    return Dict("pid" => pid, "rss_MB" => kb("VmRSS", status) / 1024, "peak_rss_MB" => kb("VmHWM", status) / 1024,
                "pss_MB" => kb("Pss", rollup) / 1024, "julia_gc_live_MB" => Base.gc_live_bytes() / 2^20)
end
