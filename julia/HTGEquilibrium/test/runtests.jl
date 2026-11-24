# Tests of the water models, HKF and the SUPCRT92 band check.
#
#   julia --project=julia/HTGEquilibrium -t 4 -e 'using Pkg; Pkg.test()'
#
# Needs the Fortran library (julia julia/HTGEquilibrium/deps/build.jl).
# Golden data: tests/reference/data/thermo.json of the Python app.
using HTGEquilibrium
using HTGEquilibrium: h2o92, jn91_epsilon, TR, PR_PA, MW_H2O_IAPWS, hkf_apparent_G
using Clapeyron
using JSON
using SHA
using TOML
using Test

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const THERMO = joinpath(REPO, "tests", "reference", "data", "thermo.json")
const GRID = joinpath(@__DIR__, "data", "pychnosz_grid.json")
const FIXTURE = joinpath(@__DIR__, "data", "hkf_ions.toml")
const H2O92D_SHA256 = "885cf390d235024db0e96c249c3e6d533cd08c31be672b26c50b90786cbc3e46"

# Python tests/helpers.py `same`/`close`: |a - b| <= max(rtol |b|, atol) (math.isclose)
# and numpy assert_allclose: |a - b| <= atol + rtol |b|.
isclose_py(a, b; rtol=1e-10, atol=1e-12) = abs(a - b) <= max(rtol * max(abs(a), abs(b)), atol)
allclose_np(a, b; rtol=1e-10, atol=1e-12) = all(abs.(a .- b) .<= atol .+ rtol .* abs.(b))
reldev(a, b) = abs(a - b) / max(abs(b), floatmin())

kelvin_pascal(T_C, P_bar) = (T_C + 273.15, P_bar * 1e5)   # as tests/test_reference_thermo.py

const RECORDS = Dict(k => HKFRecord(v) for (k, v) in TOML.parsefile(FIXTURE))
const SUPCRT = SUPCRT92Water()

@testset "HTGEquilibrium water/HKF" begin

@testset "Fortran source is the original SUPCRT92 file" begin
    @test bytes2hex(open(sha256, joinpath(@__DIR__, "..", "deps", "fortran", "H2O92D.f"))) == H2O92D_SHA256
end

@testset "SUPCRT92 golden data (thermo.json)" begin
    worst_w, worst_mu = 0.0, 0.0
    states = filter(s -> haskey(s, "electrolyte"), JSON.parsefile(THERMO; allownan=true)["states"])
    @test length(states) == 3
    for s in states
        T, P = kelvin_pascal(s["T_C"], s["P_bar"])
        w = water_properties(SUPCRT, T, P)
        for (label, ref) in s["electrolyte"]
            for (k, v) in ref["water"]
                got = getfield(w, Symbol(k))
                @test isclose_py(got, v)
                worst_w = max(worst_w, reldev(got, v))
            end
            ions = [x for x in ref["selection"] if haskey(RECORDS, x)]
            @test length(ions) == length(ref["mu0_ion"])
            mu = hkf_ion_g0(SUPCRT, [RECORDS[x] for x in ions], T, P)
            expected = Float64.(ref["mu0_ion"])
            @test allclose_np(mu, expected)
            worst_mu = max(worst_mu, maximum(reldev.(mu, expected)))
        end
    end
    @info "thermo.json: worst relative deviation" water = worst_w mu0_ion = worst_mu
end

if isfile(GRID)
    @testset "SUPCRT92 against pychnosz on the wide grid (generated, not committed)" begin
        grid = JSON.parsefile(GRID; allownan=true)
        records = [RECORDS[x] for x in grid["ions"]]
        worst_w, worst_mu, n = 0.0, 0.0, 0
        for p in grid["points"]
            T, P = kelvin_pascal(p["T_C"], p["P_bar"])
            w = water_properties(SUPCRT, T, P)
            for (k, v) in p["water"]
                v === nothing && continue
                got = getfield(w, Symbol(k))
                @test isclose_py(got, v)
                worst_w = max(worst_w, reldev(got, v))
            end
            if haskey(p, "mu0_ion")
                mu = hkf_ion_g0(SUPCRT, records, T, P)
                @test allclose_np(mu, Float64.(p["mu0_ion"]))
                worst_mu = max(worst_mu, maximum(reldev.(mu, Float64.(p["mu0_ion"]))))
            else
                @test_throws DomainError hkf_ion_g0(SUPCRT, records, T, P)
            end
            n += 1
        end
        @info "pychnosz grid: worst relative deviation" points = n water = worst_w mu0_ion = worst_mu
    end
else
    @info "No $GRID; skipping the wide pychnosz comparison (run test/dump_pychnosz.py --grid)."
end

@testset "SUPCRT92 calls" begin
    # Outside H2O92's validity range: NaN, as pychnosz.
    @test all(isnan, getfield.(Ref(water_properties(SUPCRT, 200.0, 240e5)), fieldnames(WaterProperties)))
    # Deterministic and independent of earlier calls (COMMON-block state).
    pts = [kelvin_pascal(t, p) for t in 100.0:7.5:600.0 for p in (225.0, 240.0, 300.0)]
    pts = vcat(pts, [kelvin_pascal(t, 240.0) for t in 378.0:0.05:380.0])
    forward = [water_properties(SUPCRT, T, P) for (T, P) in pts]
    backward = reverse([water_properties(SUPCRT, T, P) for (T, P) in reverse(pts)])
    @test forward == backward
    # Thread safety: the lock serialises the Fortran; results equal the serial ones.
    threaded = Vector{WaterProperties}(undef, length(pts))
    for _ in 1:5
        Threads.@threads for i in eachindex(pts)
            threaded[i] = water_properties(SUPCRT, pts[i]...)
        end
        @test threaded == forward
    end
    Threads.nthreads() == 1 && @info "Run with -t 4 to exercise the threaded test on several threads."
end

@testset "SUPCRT92 failure bands" begin
    t = @elapsed bands = supcrt92_failure_bands(240.0)
    @info "Failure bands at 240 bar ($(round(t; digits=2)) s)" bands
    # The band inside the app's electrolyte range (spike scan: 378.55-379.00 °C).
    b = only(filter(b -> b.rho == 420.04, bands))
    @test 378.50 < b.T_lo_C <= 378.55
    @test 379.00 <= b.T_hi_C < 379.05
    b2 = only(filter(b -> b.rho == 199.775, bands))
    @test 386.40 < b2.T_lo_C <= 386.45 && 386.90 <= b2.T_hi_C < 386.95
    # All bands of the spike's 0.05 °C scan are found, edges consistent with it.
    spike = Dict(225.0 => [(374.60, 374.80), (378.05, 378.40)], 240.0 => [(378.55, 379.00), (386.45, 386.90)],
                 250.0 => [(381.15, 381.75), (392.10, 392.60)], 275.0 => [(387.75, 388.45), (406.55, 406.85)],
                 300.0 => [(394.30, 395.05)])
    for (P, list) in spike, (lo, hi) in list
        @test any(b -> lo - 0.05 < b.T_lo_C <= lo && hi <= b.T_hi_C < hi + 0.05, supcrt92_failure_bands(P))
    end
    # GUI check.
    Ts = [200.0, 377.0, 378.2, 378.7, 379.3, 380.0, 386.6]
    r = supcrt92_band_check(Ts, 240.0)
    @test [x.status for x in r] == [:ok, :ok, :near_band, :in_band, :near_band, :ok, :in_band]
    @test r[4].band == b && r[3].band == b && r[1].band === nothing
    @test [x.status for x in supcrt92_band_check([378.2, 378.2], [240.0, 300.0])] == [:near_band, :ok]
    @test supcrt92_band_check([378.2], 240.0; margin_C=0.1)[1].status == :ok
    # Along the app's 240 bar sweep range, the point test fails only inside bands.
    for T in 100.0:0.25:600.0
        inside = any(b -> b.T_lo_C <= T <= b.T_hi_C, bands)
        @test supcrt92_failed(T, 240.0) == inside
    end
    # The thermo.json states are clear of all bands.
    @test all(x -> x.status == :ok, supcrt92_band_check([200.0, 300.0], 240.0))
    @test supcrt92_band_check([350.0], 250.0)[1].status == :ok
end

@testset "IAPWS-95 self-consistency" begin
    m = IAPWS95Water()
    states = [(298.15, 1e5), (473.15, 240e5), (573.15, 240e5), (623.15, 250e5), (723.15, 240e5)]
    for (T, P) in states
        w = water_properties(m, T, P)
        # Density is Clapeyron's, and is a root of p(V, T) = P.
        @test w.rho == MW_H2O_IAPWS / Clapeyron.volume(m.eos, P, T)
        @test isapprox(Clapeyron.pressure(m.eos, MW_H2O_IAPWS / w.rho, T), P; rtol=1e-9)
        @test w.epsilon == jn91_epsilon(T, w.rho / 1000)
        # Derivatives against central finite differences of the model itself.
        hT, hP = 1e-2, 1e-2 * 1e5          # K, Pa (0.01 bar)
        f(T, P) = water_properties(m, T, P)
        up, dn, pu, pd = f(T + hT, P), f(T - hT, P), f(T, P + hP), f(T, P - hP)
        @test isapprox(w.alpha, -(log(up.rho) - log(dn.rho)) / 2hT; rtol=1e-6)
        @test isapprox(w.beta, (log(pu.rho) - log(pd.rho)) / (2hP / 1e5); rtol=1e-6)
        @test isapprox(w.daldT, (up.alpha - dn.alpha) / 2hT; rtol=1e-5)
        @test isapprox(w.YBorn, (-1 / up.epsilon + 1 / dn.epsilon) / 2hT; rtol=1e-6)
        @test isapprox(w.QBorn, (-1 / pu.epsilon + 1 / pd.epsilon) / (2hP / 1e5); rtol=1e-6)
        @test isapprox(w.XBorn, (up.YBorn - dn.YBorn) / 2hT; rtol=1e-5)
    end
    # HKF runs on IAPWS-95 water and uses that model for both states.
    recs = [RECORDS[k] for k in sort(collect(keys(RECORDS)))]
    T, P = 573.15, 240e5
    mu = hkf_ion_g0(m, recs, T, P)
    @test all(isfinite, mu)
    w, wr = water_properties(m, T, P), water_properties(m, TR, PR_PA)
    @test mu == [hkf_apparent_G(r, T, P / 1e5, w, wr) * 4.184 - (r.Gf - r.Hf + 298.15 * r.Sr) * 4.184 for r in recs]
    # Thread use (Clapeyron model shared between tasks).
    out = Vector{Vector{Float64}}(undef, 8)
    Threads.@threads for i in 1:8
        out[i] = hkf_ion_g0(m, recs, T, P)
    end
    @test all(==(mu), out)
end

end
