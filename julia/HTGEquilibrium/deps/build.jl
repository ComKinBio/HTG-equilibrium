# Build the SUPCRT92 water library (H2O92) from the original Fortran source.
#
#   julia julia/HTGEquilibrium/deps/build.jl
#
# (`Pkg.build("HTGEquilibrium")` runs this file too.) It compiles
# vendor/supcrt92/H2O92D.f, which is the unmodified SUPCRT92 file, into
# deps/build/libh2o92.<so|dylib|dll>. The output is build output and is not
# committed. Set FC to use a Fortran compiler other than `gfortran`.
#
# Flags:
#   -O2                 the optimisation level pychnosz's library gives results
#                       identical to (checked bit for bit on the test grid)
#   -ffp-contract=off   no fused multiply-add, so results do not depend on the
#                       CPU (aarch64 compilers contract by default)
#   -std=legacy         the file is Fortran 77 with legacy constructs; this
#                       only silences the warnings about them
# Cross-platform builds (Windows, macOS, a JLL through BinaryBuilder) are an
# M5 packaging task.

using Libdl

const SRC = joinpath(@__DIR__, "..", "vendor", "supcrt92", "H2O92D.f")
const OUTDIR = joinpath(@__DIR__, "build")
const LIB = joinpath(OUTDIR, "libh2o92." * Libdl.dlext)

function build_h2o92(; fc=get(ENV, "FC", "gfortran"))
    compiler = Sys.which(fc)
    compiler === nothing && error("Fortran compiler `$fc` not found. Install gfortran " *
                                  "(e.g. `apt install gfortran`) or set FC.")
    mkpath(OUTDIR)
    cmd = `$compiler -O2 -fPIC -shared -ffp-contract=off -std=legacy -o $LIB $SRC`
    @info "Building SUPCRT92 H2O92" cmd
    run(cmd)
    return LIB
end

build_h2o92()
