#!/bin/bash
# Start-up time and idle memory, N runs each without and with the sysimage.
#   bash measure.sh [N]
# Each run: launch -> window usable -> 10 s idle -> memory -> window closed -> check no
# julia/electron process of the run is left.
cd "$(dirname "$0")"
N=${1:-3}
run() {
  local label=$1; shift
  for i in $(seq 1 "$N"); do
    load=$(cut -d' ' -f1 /proc/loadavg)
    out=$(M0_T0=$(date +%s.%N) M0_IDLE=10 M0_AUTOCLOSE=1 julia --project=. "$@" run.jl 2>&1)
    usable=$(grep -o '[0-9.]* s  USABLE' <<<"$out" | cut -d' ' -f1)
    page=$(grep -o '[0-9.]* s  window shows' <<<"$out" | cut -d' ' -f1)
    mem=$(grep -o 'memory: .*' <<<"$out")
    sleep 1
    left=$(pgrep -f "m0_window/run.jl|jlel-" | wc -l)
    echo "$label run $i (load $load): page loaded ${page} s, usable ${usable} s; $mem; processes left: $left"
  done
}
run "no sysimage"
[ -f m0_sysimage.so ] && run "sysimage" -J m0_sysimage.so
