# Tests of HTGEquilibrium. One file per area; each compares with the Python app.
#
#   julia --project=julia/HTGEquilibrium -t 4 -e 'using Pkg; Pkg.test()'
#
# The water tests need the Fortran library (julia julia/HTGEquilibrium/deps/build.jl).
using Test

@testset "HTGEquilibrium" begin
    include("test_formulas.jl")
    include("test_water.jl")
    include("test_units.jl")
    include("test_number_lists.jl")
    include("test_settings.jl")
    include("test_case.jl")
    include("test_feed.jl")
    include("test_databases.jl")
    include("test_registry.jl")
    include("test_standard_states.jl")
    include("test_phases.jl")
    include("test_optimizer.jl")
    include("test_equilibrium.jl")
    include("test_results.jl")
    include("test_sweep.jl")
    include("test_queue.jl")
    include("test_cli.jl")
end
