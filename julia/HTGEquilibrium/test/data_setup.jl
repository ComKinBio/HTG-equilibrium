# Shared setup of the data-layer tests (test_databases.jl, test_registry.jl,
# test_standard_states.jl): paths, the two reference models and Python's
# comparison rules. Included once into Main; the test files use it as `Main.DataTestSetup`.
module DataTestSetup

using HTGEquilibrium
using HTGEquilibrium: _exception_text
using JSON
using HTGEquilibrium: OrderedDict   # OrderedCollections is not a test dependency

export REPO, FROZEN, SHIPPED, REFERENCE_DATA, MODELS, same, isclose_py, allclose_np, reldev, plain,
       outcome, reference

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
# The databases and species.toml the reference results were generated with
# (tests/reference/README.md); the shipped data/ folder may change.
const FROZEN = joinpath(REPO, "tests", "reference", "databases")
const SHIPPED = joinpath(REPO, "data")
const REFERENCE_DATA = joinpath(REPO, "tests", "reference", "data")

reference(name) = JSON.parsefile(joinpath(REFERENCE_DATA, "$name.json"); allownan=true)

# As tests/conftest.py: 'base' without extra NASA names, 'acetone' with the species.toml.
const MODELS = Dict("base" => load_model(FROZEN; nasa_names=Dict()), "acetone" => load_model(FROZEN))

# Python tests/helpers.py: `same` (math.isclose; NaN equals NaN; None equals NaN) and
# `close` (numpy assert_allclose: |a - b| <= atol + rtol |b|).
isclose_py(a, b; rtol=1e-10, atol=1e-12) = abs(a - b) <= max(rtol * max(abs(a), abs(b)), atol)
allclose_np(a, b; rtol=1e-10, atol=1e-12) =
    length(a) == length(b) && all(abs.(a .- b) .<= atol .+ rtol .* abs.(b))
reldev(a, b) = a == b ? 0.0 : abs(a - b) / max(abs(b), floatmin())

_nan(x) = x === nothing || x === missing || (x isa AbstractFloat && isnan(x))
function same(a, b; rtol=1e-10, atol=1e-12)
    if a isa AbstractDict && b isa AbstractDict
        return Set(String.(keys(a))) == Set(String.(keys(b))) &&
               all(same(a[k], b[k]; rtol, atol) for k in keys(a))
    elseif (a isa AbstractVector || a isa Tuple) && (b isa AbstractVector || b isa Tuple)
        return length(a) == length(b) && all(same(x, y; rtol, atol) for (x, y) in zip(a, b))
    elseif _nan(a) && _nan(b)
        return true
    elseif a isa Bool || b isa Bool
        return (a isa Real && b isa Real) && Bool(a) == Bool(b)
    elseif a isa Real && b isa Real
        return isclose_py(Float64(a), Float64(b); rtol, atol)
    end
    return isequal(a, b)
end

"""Julia values as the Python records they stand for (rows as dicts, tuples as lists)."""
plain(r::CatalogueRow) = OrderedDict{String,Any}(String(k) => getfield(r, k) for k in fieldnames(CatalogueRow))
plain(t::NamedTuple) = OrderedDict{String,Any}(String(k) => plain(v) for (k, v) in pairs(t))
plain(d::AbstractDict) = OrderedDict{String,Any}(String(k) => plain(v) for (k, v) in d)
plain(v::Union{AbstractVector,Tuple}) = Any[plain(x) for x in v]
plain(x) = x

"""As tests/test_reference_thermo.py `outcome`: the value, or 'ExceptionType: message' (Python names)."""
function outcome(f)
    try
        return f()
    catch err
        err isa ArgumentError && return "ValueError: " * err.msg
        err isa KeyError && return "KeyError: " * _exception_text(err)
        rethrow()
    end
end

end # module
