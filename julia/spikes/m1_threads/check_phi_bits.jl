# Why 5 of 53 fugacity-coefficient vectors differ from Python's in the last bits.
# Hypothesis: juliacall passes numpy arrays as PythonCall.PyArray, an AbstractArray
# that is not a DenseArray, so Clapeyron's reductions (sum, dot) take Julia's
# generic (non-SIMD) path; the replay passes a Vector{Float64}. A plain wrapper
# type that is only an AbstractVector reproduces that path.
#   julia --project=. check_phi_bits.jl
include(joinpath(@__DIR__, "replay.jl"))
using .Replay

struct Plain{T} <: AbstractVector{T}
    v::Vector{T}
end
Base.size(p::Plain) = size(p.v)
Base.getindex(p::Plain, i::Int) = p.v[i]
Base.IndexStyle(::Type{<:Plain}) = IndexLinear()

cases = Replay.load_cases(joinpath(@__DIR__, "problems.json"))
let nvec = 0, vec_eq = 0, plain_eq = 0
    for c in cases, (name, ph) in Replay.build_phases(c)
        rec = c["phases"][name]
        haskey(rec, "check_n") || continue
        ref = Float64.(rec["check_phi"])
        n = Float64.(rec["check_n"])
        if ph isa Replay.FluidPhase
            a = Replay.fugacity_coefficient(ph.model, ph.P, ph.T, n; phase = ph.root)
            b = Replay.fugacity_coefficient(ph.model, ph.P, ph.T, Plain(n); phase = ph.root)
        else  # phi() perturbs n before the call; reproduce its input
            a = Replay.phi(ph, n)
            m = copy(n)  # same probe arithmetic, then call with the wrapper
            s, z = ph.split, ph.z
            imb = Replay.dot(m[s+1:end], z)
            sm = hypot(imb, 2e-14)
            p1, n1 = findfirst(>(0), z), findfirst(<(0), z)
            m[s+p1] += (sm - imb) / (2 * z[p1]); m[s+n1] += (sm + imb) / (-2 * z[n1])
            m[1:s] .+= 1e-14; m[s+1:end] .+= ph.balanced .* 1e-14
            b = Replay.fugacity_coefficient(ph.model, ph.P, ph.T, Plain(m); phase = "liquid")
        end
        nvec += 1; vec_eq += a == ref; plain_eq += b == ref
        a == ref || println(c["name"], " ", name, ": Vector differs, AbstractVector wrapper ",
                            b == ref ? "matches Python bitwise" : "also differs")
    end
    println("bitwise equal to Python: Vector input $vec_eq/$nvec, non-dense AbstractVector input $plain_eq/$nvec")
end
