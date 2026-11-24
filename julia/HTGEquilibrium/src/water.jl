# Water models: the properties of pure water that the electrolyte phase and
# HKF need. Python counterpart: `standard_states.supcrt_water`.
#
#   SUPCRT92Water()  the original SUPCRT92 Fortran (H2O92: Haar-Gallagher-Kell
#                    and Levelt Sengers equations of state, Johnson & Norton
#                    1991 dielectric constant), called with ccall. This is what
#                    pychnosz runs, so it reproduces the Python app exactly.
#   IAPWS95Water()   IAPWS-95 density from Clapeyron, with the JN91 dielectric
#                    constant and Born functions computed here. A different
#                    model: results differ from SUPCRT92 (see the HKF spike).

# Reference state of HKF and SUPCRT92.
const TR = 298.15          # K
const PR = 1.0             # bar
const PR_PA = 1.0e5        # Pa

"""
Water properties at one (T, P), with the fields and units of the Python
`standard_states.supcrt_water`:

| field   | unit    | meaning                                  |
|---------|---------|------------------------------------------|
| rho     | kg/m3   | density                                  |
| epsilon | –       | relative permittivity (dielectric const.)|
| QBorn   | 1/bar   | Born Q = (d(-1/epsilon)/dP)_T            |
| XBorn   | 1/K^2   | Born X = (dY/dT)_P                       |
| YBorn   | 1/K     | Born Y = (d(-1/epsilon)/dT)_P            |
| alpha   | 1/K     | isobaric expansivity, -(d ln rho/dT)_P   |
| daldT   | 1/K^2   | (d alpha/dT)_P                           |
| beta    | 1/bar   | isothermal compressibility, (d ln rho/dP)_T |
"""
struct WaterProperties
    rho::Float64
    epsilon::Float64
    QBorn::Float64
    XBorn::Float64
    YBorn::Float64
    alpha::Float64
    daldT::Float64
    beta::Float64
end

const WATER_FIELDS = fieldnames(WaterProperties)

"""Water properties as a `Dict{String,Float64}` with the Python key names."""
water_dict(w::WaterProperties) = Dict(String(k) => getfield(w, k) for k in WATER_FIELDS)

const NAN_WATER = WaterProperties(ntuple(_ -> NaN, length(WATER_FIELDS))...)

"""Abstract type of the selectable water models."""
abstract type WaterModel end

"""
    water_properties(model::WaterModel, T, P) -> WaterProperties

Water properties at temperature `T` [K] and pressure `P` [Pa] (the units of
the Python `supcrt_water`).
"""
function water_properties end

# ---------------------------------------------------------------------------
# SUPCRT92: the original Fortran (vendor/supcrt92/H2O92D.f)
# ---------------------------------------------------------------------------

"""
    SUPCRT92Water()

SUPCRT92 water (H2O92D.f, J. W. Johnson 1991) through its original Fortran,
compiled by `deps/build.jl`. Same settings as pychnosz/CHNOSZ: T in °C, P in
bar, density in g/cm3, Levelt Sengers equation near the critical point, JN91
dielectric constant.

Thread safety: H2O92 keeps its state in COMMON blocks and SAVEd local
variables, so it is not reentrant. Every call goes through one global lock
(`H2O92_LOCK`); calls from several threads are therefore safe but run one at
a time.

Outside the validity range of H2O92 the result is all NaN, as in pychnosz.
"""
struct SUPCRT92Water <: WaterModel end

"""Lock that serialises all calls into the (non-reentrant) H2O92 Fortran."""
const H2O92_LOCK = ReentrantLock()
const _H2O92_FN = Ref{Ptr{Cvoid}}(C_NULL)
# LVS region limits read from the library when it is loaded.
const _LVS_TLIMITS = Ref((NaN, NaN))   # Tmin1, Tmax [K]
const _LVS_DLIMITS = Ref((NaN, NaN))   # Dmin, Dmax [kg/m3]

"""Path of the compiled H2O92 library (`HTG_H2O92_LIB` overrides the default)."""
h2o92_library_path() = get(ENV, "HTG_H2O92_LIB",
                           joinpath(dirname(@__DIR__), "deps", "build", "libh2o92." * Libdl.dlext))

# Call with H2O92_LOCK held.
function _h2o92_fn()
    if _H2O92_FN[] == C_NULL
        path = h2o92_library_path()
        isfile(path) || error("SUPCRT92 library not found at $path. Build it with " *
                              "`julia $(joinpath(dirname(@__DIR__), "deps", "build.jl"))` (needs gfortran).")
        handle = Libdl.dlopen(path)
        # COMMON /coefs/ a(20), q(20), x(11): x(1) = Tmin1 and x(3) = Tmax [K],
        # x(4) = Dmin and x(5) = Dmax [kg/m3] of the LVS (critical) region.
        coefs = Ptr{Float64}(Libdl.dlsym(handle, :coefs_))
        _LVS_TLIMITS[] = (unsafe_load(coefs, 41), unsafe_load(coefs, 43))
        _LVS_DLIMITS[] = (unsafe_load(coefs, 44), unsafe_load(coefs, 45))
        _H2O92_FN[] = Libdl.dlsym(handle, :h2o92_)
    end
    return _H2O92_FN[]
end

"""
    h2o92(T_C, x; iopt=2) -> NamedTuple

One raw call of SUBROUTINE H2O92 for a single phase (`isat = 0`), with the
specs pychnosz/CHNOSZ use: `(it, id, ip, ih, itripl, isat, iopt, useLVS,
epseqn, icrit) = (2, 2, 2, 5, 1, 0, iopt, 1, 4, 0)`.
`iopt = 2`: `x` is the pressure in bar; `iopt = 1`: `x` is the density in g/cm3.
Returns `(T_C, P_bar, D, props, error, lvs)`: the state H2O92 returns, the 46
property values (vapour/liquid interleaved; for a single phase use the odd
indices), the error flag and whether the Levelt Sengers (critical region)
equation was used.
"""
function h2o92(T_C::Real, x::Real; iopt::Integer=2)
    specs = Int32[2, 2, 2, 5, 1, 0, iopt, 1, 4, 0]
    states = iopt == 2 ? Float64[T_C, x, 1.0, 0.0] : Float64[T_C, 0.0, x, 0.0]
    props = zeros(Float64, 46)
    err = Ref{Int32}(0)   # Fortran default LOGICAL is 4 bytes
    @lock H2O92_LOCK begin
        fn = _h2o92_fn()
        ccall(fn, Cvoid, (Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ref{Int32}),
              specs, states, props, err)
    end
    return (T_C=states[1], P_bar=states[2], D=states[3], props=props,
            error=err[] != 0, lvs=specs[10] == 1)
end

# Index of each property in H2O92's props(46) for the single-phase (first) state.
# Order in H2O92: A G S U H Cv Cp Speed alpha beta diel visc tcond surten tdiff
# Prndtl visck albe ZBorn YBorn QBorn daldT XBorn.
const _H2O92_INDEX = (alpha=17, beta=19, epsilon=21, YBorn=39, QBorn=41, daldT=43, XBorn=45)

function water_properties(::SUPCRT92Water, T::Real, P::Real)
    # Same arithmetic as standard_states.supcrt_water + pychnosz.
    r = h2o92(Float64(T) - 273.15, Float64(P) / 100000.0)
    r.error && return NAN_WATER
    p, i = r.props, _H2O92_INDEX
    return WaterProperties(r.D * 1000, p[i.epsilon], p[i.QBorn], p[i.XBorn], p[i.YBorn],
                           p[i.alpha], p[i.daldT], p[i.beta])
end

# ---------------------------------------------------------------------------
# Johnson & Norton (1991) dielectric constant and Born functions (used with
# IAPWS-95). Ported from the HKF spike, which reproduces SUPCRT92's JN91 to
# 2e-16 given the same density. H2O92D.f writes the JN91 coefficients as
# single-precision literals (0.1470333593E+02), so SUPCRT92 uses their
# Float32 values; they are rounded the same way here.
# ---------------------------------------------------------------------------
const JN91_A = Float64.(Float32.((0.1470333593e2, 0.2128462733e3, -0.1154445173e3, 0.1955210915e2,
                                  -0.8330347980e2, 0.3213240048e2, -0.6694098645e1, -0.3786202045e2,
                                  0.6887359646e2, -0.2729401652e2)))

"""JN91 dielectric constant at T [K] and density D [g/cm3]."""
function jn91_epsilon(T, D)
    a = JN91_A
    th = T / 298.15
    c2 = a[1] / th
    c3 = a[2] / th + a[3] + a[4] * th
    c4 = a[5] / th + a[6] * th + a[7] * th^2
    c5 = a[8] / th^2 + a[9] / th + a[10]
    return 1 + c2 * D + c3 * D^2 + c4 * D^3 + c5 * D^4
end

"""
    water_from_density(T, rho, alpha, daldT, beta) -> WaterProperties

Dielectric constant (JN91) and Born functions from T [K], the density
[kg/m3] and its derivatives (alpha 1/K, daldT 1/K2, beta 1/bar).
"""
function water_from_density(T, rho, alpha, daldT, beta)
    D = rho / 1000
    eps = jn91_epsilon(T, D)
    g = ForwardDiff.gradient(x -> jn91_epsilon(x[1], x[2]), [T, D])
    H = ForwardDiff.hessian(x -> jn91_epsilon(x[1], x[2]), [T, D])
    eT, eD = g
    eTT, eTD, eDD = H[1, 1], H[1, 2], H[2, 2]
    dDdT = -D * alpha
    dDdP = D * beta
    d2DdT2 = D * (alpha^2 - daldT)
    dedT = eT + eD * dDdT
    dedP = eD * dDdP
    d2edT2 = eTT + 2 * eTD * dDdT + eDD * dDdT^2 + eD * d2DdT2
    Q = dedP / eps^2
    Y = dedT / eps^2
    X = d2edT2 / eps^2 - 2 * dedT^2 / eps^3
    return WaterProperties(rho, eps, Q, X, Y, alpha, daldT, beta)
end

# ---------------------------------------------------------------------------
# IAPWS-95 (Clapeyron)
# ---------------------------------------------------------------------------

"""
    IAPWS95Water()

IAPWS-95 water from Clapeyron (`Clapeyron.IAPWS95()`), with the JN91
dielectric constant and Born functions. alpha, beta and dalpha/dT are exact
derivatives (ForwardDiff through two Newton steps on Clapeyron's
`pressure(V, T)` from the converged volume).

This is a different equation of state from SUPCRT92 (HGK/LVS): density differs
by 1e-4 to 6e-3 relative where the app uses HKF, so it does not reproduce the
Python golden data. The 25 °C / 1 bar reference values used by HKF are also
taken from IAPWS-95.
"""
struct IAPWS95Water{M} <: WaterModel
    eos::M
end
IAPWS95Water() = IAPWS95Water(Clapeyron.IAPWS95())

const MW_H2O_IAPWS = 18.015268e-3   # kg/mol

# Two Newton steps on p(V, T) = P from the converged V0. With Dual T and P the
# result carries exact first and second derivatives of the root.
function _newton_volume(eos, V0, T, P_Pa)
    V = V0
    for _ in 1:2
        f = Clapeyron.pressure(eos, V, T) - P_Pa
        dfdV = ForwardDiff.derivative(v -> Clapeyron.pressure(eos, v, T), V)
        V = V - f / dfdV
    end
    return V
end

function water_properties(m::IAPWS95Water, T::Real, P::Real)
    T, P = Float64(T), Float64(P)
    V0 = Clapeyron.volume(m.eos, P, T)
    lnrho(x) = log(MW_H2O_IAPWS / _newton_volume(m.eos, V0, x[1], x[2] * 1e5))
    x = [T, P / 1e5]
    g = ForwardDiff.gradient(lnrho, x)
    H = ForwardDiff.hessian(lnrho, x)
    return water_from_density(T, MW_H2O_IAPWS / V0, -g[1], -H[1, 1], g[2])
end
