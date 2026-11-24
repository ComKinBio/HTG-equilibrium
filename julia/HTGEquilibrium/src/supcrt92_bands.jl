# SUPCRT92 failure bands near the pseudo-critical line, for the GUI warning.
#
# H2O92 solves for the density at (T, P) with the Levelt Sengers (LVS)
# equation inside the critical region (643-695 K). Its solver, denLVS, takes at
# most 20 Newton steps and clamps the density to the LVS limits Dmin = 199.775
# and Dmax = 420.04 kg/m3. In bands about 0.5 °C wide on each side of the
# pseudo-critical line it does not converge and returns one of the clamp
# values, so the state is up to a few bar off the requested pressure. A few
# isolated points elsewhere in the LVS region stop after 20 steps without
# reaching a clamp value. H2O92 reports neither as an error. The result is
# deterministic, so the golden data (and pychnosz) contain these values.
#
# Detection (no table, works at any pressure):
#   1. Call H2O92 at (T, P) and take the density D it returns.
#   2. Call H2O92 again at (T, D) (iopt = 1) to get the pressure P(T, D) of
#      that density on the same equation of state.
#   3. The point has failed if D is exactly Dmin or Dmax, or if the density
#      error implied by the pressure mismatch, beta * |P(T, D) - P|, is above
#      1e-6 (relative).
# Converged points have an implied density error below 1e-6 (99 % below
# 1e-10); see README for the scan behind this threshold.
#
# Bands at one pressure are found by applying this test every 0.01 °C over the
# temperature range of the LVS region (read from the library: Tmin1 to Tmax of
# COMMON /coefs/) and refining each edge by bisection to 1e-4 °C. A scan of
# 100-600 °C x 150-350 bar every 0.05 °C (5 bar) found no failures outside the
# LVS region. `supcrt92_band_check` also tests every requested temperature
# directly, so a failing temperature is reported even if it lies outside the
# scanned range or in a band narrower than the scan step.

const FAILURE_DLNRHO = 1e-6
const BAND_SCAN_STEP_C = 0.01
const BAND_EDGE_TOL_C = 1e-4

"""
A temperature interval at one pressure where SUPCRT92 (H2O92) does not
converge. `T_lo_C` and `T_hi_C` are the first and last failing temperatures
(to 1e-4 °C); `rho` is the density H2O92 returns in the middle of the band
[kg/m3]: 420.04 or 199.775 for the clamp bands.
"""
struct SUPCRT92Band
    P_bar::Float64
    T_lo_C::Float64
    T_hi_C::Float64
    rho::Float64
end

"""
    supcrt92_failed(T_C, P_bar) -> Bool

Whether SUPCRT92's density solver fails at `T_C` [°C], `P_bar` [bar]: it
returns an LVS clamp density (199.775 or 420.04 kg/m3) or a density whose
pressure on the same equation of state is off by more than a relative density
error of 1e-6. States outside H2O92's validity range return `false`.
"""
function supcrt92_failed(T_C::Real, P_bar::Real)
    r = h2o92(Float64(T_C), Float64(P_bar))
    r.error && return false
    Dmin, Dmax = _LVS_DLIMITS[]
    (r.D == Dmin / 1000 || r.D == Dmax / 1000) && return true
    back = h2o92(Float64(T_C), r.D; iopt=1)
    beta = r.props[_H2O92_INDEX.beta]
    return !(beta * abs(back.P_bar - P_bar) <= FAILURE_DLNRHO)
end

const _BANDS_LOCK = ReentrantLock()
const _BANDS_CACHE = Dict{Float64,Vector{SUPCRT92Band}}()

# LVS temperature range [°C] (COMMON /coefs/ x(1) = Tmin1, x(3) = Tmax).
function _lvs_range_C()
    @lock H2O92_LOCK _h2o92_fn()
    Tmin, Tmax = _LVS_TLIMITS[]
    return Tmin - 273.15, Tmax - 273.15
end

function _edge(ok_T, fail_T, P_bar)
    while abs(fail_T - ok_T) > BAND_EDGE_TOL_C
        mid = (ok_T + fail_T) / 2
        supcrt92_failed(mid, P_bar) ? (fail_T = mid) : (ok_T = mid)
    end
    return fail_T
end

function _scan_bands(P_bar::Float64)
    lo, hi = _lvs_range_C()
    Ts = collect(range(lo, hi; step=BAND_SCAN_STEP_C))
    failed = [supcrt92_failed(T, P_bar) for T in Ts]
    bands = SUPCRT92Band[]
    i = 1
    while i <= length(Ts)
        if !failed[i]
            i += 1
            continue
        end
        j = i
        while j < length(Ts) && failed[j + 1]
            j += 1
        end
        T_lo = i > 1 ? _edge(Ts[i - 1], Ts[i], P_bar) : Ts[i]
        T_hi = j < length(Ts) ? _edge(Ts[j + 1], Ts[j], P_bar) : Ts[j]
        mid = h2o92((T_lo + T_hi) / 2, P_bar)
        push!(bands, SUPCRT92Band(P_bar, T_lo, T_hi, mid.D * 1000))
        i = j + 1
    end
    return bands
end

"""
    supcrt92_failure_bands(P_bar) -> Vector{SUPCRT92Band}

All SUPCRT92 failure bands at `P_bar` [bar] in the LVS temperature range
(369.85-421.85 °C). Computed on first use for each pressure (~0.3 s) and
cached. Safe to call from several threads.
"""
function supcrt92_failure_bands(P_bar::Real)
    P = Float64(P_bar)
    cached = @lock _BANDS_LOCK get(_BANDS_CACHE, P, nothing)
    cached === nothing || return cached
    bands = _scan_bands(P)
    @lock _BANDS_LOCK _BANDS_CACHE[P] = bands
    return bands
end

"""
    supcrt92_band_check(T_C, P_bar; margin_C=0.5)

For the GUI warning before a sweep with the SUPCRT92 water model. `T_C` is a
vector of temperatures [°C]; `P_bar` is one pressure [bar] or a vector with one
pressure per temperature. Returns one named tuple per temperature:

    (T_C, P_bar, status, band)

with `status`
- `:in_band`   SUPCRT92 fails at this temperature (inside a band, or the point
               test fails directly);
- `:near_band` within `margin_C` of a band;
- `:ok`        otherwise,

and `band` the nearest `SUPCRT92Band` for `:in_band`/`:near_band` (`nothing`
if the point failed outside any scanned band or the status is `:ok`).
"""
function supcrt92_band_check(T_C::AbstractVector{<:Real}, P_bar; margin_C::Real=0.5)
    Ps = P_bar isa Real ? fill(Float64(P_bar), length(T_C)) : Float64.(collect(P_bar))
    length(Ps) == length(T_C) || throw(DimensionMismatch("need one pressure, or one per temperature"))
    return map(zip(Float64.(T_C), Ps)) do (T, P)
        bands = supcrt92_failure_bands(P)
        nearest, dist = nothing, Inf
        for b in bands
            d = T < b.T_lo_C ? b.T_lo_C - T : T > b.T_hi_C ? T - b.T_hi_C : 0.0
            d < dist && ((nearest, dist) = (b, d))
        end
        status = (dist == 0 || supcrt92_failed(T, P)) ? :in_band :
                 dist <= margin_C ? :near_band : :ok
        (T_C=T, P_bar=P, status=status, band=(status != :ok && dist <= margin_C) ? nearest : nothing)
    end
end
