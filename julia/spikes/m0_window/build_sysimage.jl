# Build m0_sysimage.so (not committed). Needs a desktop session: the precompile
# workload opens the Electron window.
#   julia --project=. build_sysimage.jl
using PackageCompiler

create_sysimage([:M0App, :Genie, :Stipple, :StippleUI, :StipplePlotly, :Electron, :DataFrames, :JSON, :HTTP];
                sysimage_path = joinpath(@__DIR__, "m0_sysimage.so"),
                precompile_execution_file = joinpath(@__DIR__, "precompile.jl"))
