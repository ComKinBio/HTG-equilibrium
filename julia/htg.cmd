@echo off
rem htg: the HTG equilibrium command line, Julia version (Windows).
rem
rem     julia\htg.cmd check cases\manure_gas_min.toml
rem     julia\htg.cmd run cases\manure_gas_min.toml --workers 4
rem     julia\htg.cmd gui
rem
rem Runs "julia -m HTGEquilibrium ARGS..." with this repository's package
rem (julia\HTGEquilibrium) and one thread per CPU thread, so that --workers N
rem can use N threads. The command gui runs "julia -m HTGApp ARGS..." with the
rem app package (julia\HTGApp) instead, so that the other commands do not load
rem the GUI packages. JULIA_NUM_THREADS chooses another thread count, JULIA
rem another julia executable (default: julia on PATH).
setlocal
set "HTG_THREADS=%JULIA_NUM_THREADS%"
if "%HTG_THREADS%"=="" set "HTG_THREADS=auto"
set "HTG_JULIA=%JULIA%"
if "%HTG_JULIA%"=="" set "HTG_JULIA=julia"
set "HTG_PACKAGE=HTGEquilibrium"
call :find_command %*
"%HTG_JULIA%" --project="%~dp0%HTG_PACKAGE%" --threads=%HTG_THREADS% --startup-file=no -m %HTG_PACKAGE% %*
exit /b %ERRORLEVEL%

rem The command is the first word that is not an option of htg itself (-h, --data DIR).
:find_command
if "%~1"=="" goto :eof
set "HTG_ARG=%~1"
if "%HTG_ARG%"=="--d" goto skip_value
if "%HTG_ARG%"=="--da" goto skip_value
if "%HTG_ARG%"=="--dat" goto skip_value
if "%HTG_ARG%"=="--data" goto skip_value
if "%HTG_ARG:~0,1%"=="-" goto next_word
if "%HTG_ARG%"=="gui" set "HTG_PACKAGE=HTGApp"
goto :eof
:skip_value
shift
:next_word
shift
goto find_command
