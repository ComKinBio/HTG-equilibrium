# M0 spike launcher: Genie server on a free local port + Electron.jl window.
#
#   julia --project=. run.jl                       # window stays until you close it
#   julia --project=. -J m0_sysimage.so run.jl     # same with the sysimage (build_sysimage.jl)
#
# Environment variables:
#   M0_T0          launch time as Unix seconds (date +%s.%N), so Julia's own start-up counts too
#   M0_SELFTEST=1  drive every component in the window and print PASS/FAIL per component
#   M0_IDLE=s      sit idle s seconds after start-up, then print the memory of all processes
#   M0_AUTOCLOSE=1 close the window (as the user would) at the end
#   M0_SHOT=path   write a PNG screenshot of the window after start-up
#   M0_GTK_DIALOGS=1  use Electron's GTK file dialogs instead of the XDG portal (GNOME)
#   M0_JSDIR=dir   debugging: run each *.js dropped into dir in the page (*.main.js in Electron's
#                  main process) and write the result to *.out
t0 = parse(Float64, get(ENV, "M0_T0", string(time())))
println("[M0] $(lpad(round(time() - t0; digits = 2), 6)) s  julia running"); flush(stdout)

ENV["GENIE_ENV"] = get(ENV, "GENIE_ENV", "prod")   # minified Vue/Quasar, no file watching
using M0App
M0App.T0[] = t0
M0App.stamp("packages loaded")
M0App.main(t0)
exit(0)
