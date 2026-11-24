@echo off
rem htg: the HTG equilibrium command line, Julia version (Windows).
rem
rem     julia\htg.cmd check cases\manure_gas_min.toml
rem     julia\htg.cmd run cases\manure_gas_min.toml --workers 4
rem
rem Runs "julia -m HTGEquilibrium ARGS..." with this repository's package
rem (julia\HTGEquilibrium) and one thread per CPU thread, so that --workers N
rem can use N threads. JULIA_NUM_THREADS chooses another thread count, JULIA
rem another julia executable (default: julia on PATH).
setlocal
set "HTG_THREADS=%JULIA_NUM_THREADS%"
if "%HTG_THREADS%"=="" set "HTG_THREADS=auto"
set "HTG_JULIA=%JULIA%"
if "%HTG_JULIA%"=="" set "HTG_JULIA=julia"
"%HTG_JULIA%" --project="%~dp0HTGEquilibrium" --threads=%HTG_THREADS% --startup-file=no -m HTGEquilibrium %*
exit /b %ERRORLEVEL%
