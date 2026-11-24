# Option 2 of the Phase 2 plan (section 4): IAPWS-95 water density from
# Clapeyron, plus the parts of SUPCRT92 the app needs, written in Julia:
#   - Johnson & Norton (1991) dielectric constant (SUPCRT92 H2O92D.f, JN91)
#   - Born functions Q, Y, X (derivatives of Z = -1/epsilon)
#   - g-function of Shock et al. (1992) for omega(T, P) of charged species
#   - revised HKF standard Gibbs energy (Tanger & Helgeson 1988; Shock et al. 1992)
# Units follow SUPCRT92/CHNOSZ: T in K, P in bar, density in g/cm3 inside the
# formulas, energies in cal/mol (the app multiplies by 4.184).
module HKFWater

using Clapeyron
using ForwardDiff

export WaterState, iapws95_water, water_from_density, hkf_G, hkf_records_from_dict, HKFRecord

const MW_H2O = 18.015268e-3        # kg/mol (IAPWS-95)
const TR = 298.15
const PR = 1.0
const THETA = 228.0
const PSI = 2600.0
const ETA = 1.66027e5              # Angstrom cal/mol

const _MODEL = Ref{Any}(nothing)
iapws() = (_MODEL[] === nothing && (_MODEL[] = IAPWS95()); _MODEL[])

# ---------------------------------------------------------------------------
# Johnson & Norton (1991) dielectric constant, as coded in SUPCRT92 (JN91).
# D: density in g/cm3. The DATA statement in H2O92D.f writes the coefficients
# as single-precision literals (E+02, not D+02), so SUPCRT92 uses their
# Float32 values; rounding them the same way reproduces pychnosz's epsilon to
# 2e-16 (with the full decimal values it is off by ~1e-7).
# ---------------------------------------------------------------------------
const JN_A = Float64.(Float32.((0.1470333593e2, 0.2128462733e3, -0.1154445173e3, 0.1955210915e2,
              -0.8330347980e2, 0.3213240048e2, -0.6694098645e1, -0.3786202045e2,
              0.6887359646e2, -0.2729401652e2)))

function jn91_epsilon(T, D)
    a = JN_A
    th = T / 298.15
    c1 = one(th)
    c2 = a[1] / th
    c3 = a[2] / th + a[3] + a[4] * th
    c4 = a[5] / th + a[6] * th + a[7] * th^2
    c5 = a[8] / th^2 + a[9] / th + a[10]
    return c1 + c2 * D + c3 * D^2 + c4 * D^3 + c5 * D^4
end

"""Water state as `standard_states.supcrt_water` returns it (rho kg/m3, alpha 1/K,
daldT 1/K2, beta 1/bar, Born Q 1/bar, Y 1/K, X 1/K2)."""
struct WaterState
    rho::Float64
    epsilon::Float64
    QBorn::Float64
    XBorn::Float64
    YBorn::Float64
    alpha::Float64
    daldT::Float64
    beta::Float64
end

"""Dielectric constant and Born functions from T and the density with its
derivatives (alpha = -(dlnD/dT)_P, beta = (dlnD/dP)_T, daldT = (dalpha/dT)_P)."""
function water_from_density(T, rho_kgm3, alpha, daldT, beta)
    D = rho_kgm3 / 1000
    eps = jn91_epsilon(T, D)
    gT = ForwardDiff.gradient(x -> jn91_epsilon(x[1], x[2]), [T, D])
    H = ForwardDiff.hessian(x -> jn91_epsilon(x[1], x[2]), [T, D])
    eT, eD = gT
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
    return WaterState(rho_kgm3, eps, Q, X, Y, alpha, daldT, beta)
end

# ---------------------------------------------------------------------------
# IAPWS-95 density (Clapeyron) with AD-exact T and P derivatives: solve once in
# Float64, then take two Newton steps in Dual numbers (first and second
# derivatives of the root are then exact).
# ---------------------------------------------------------------------------
function iapws_volume(T::Float64, P_Pa::Float64)
    return Clapeyron.volume(iapws(), P_Pa, T)
end

function _newton_volume(V0, T, P_Pa)
    model = iapws()
    V = V0
    for _ in 1:2
        f = Clapeyron.pressure(model, V, T) - P_Pa
        dfdV = ForwardDiff.derivative(v -> Clapeyron.pressure(model, v, T), V)
        V = V - f / dfdV
    end
    return V
end

"""Water properties from IAPWS-95 (Clapeyron) + JN91 at T [K], P [bar]."""
function iapws95_water(T::Float64, P_bar::Float64)
    V0 = iapws_volume(T, P_bar * 1e5)
    lnrho(x) = log(MW_H2O / _newton_volume(V0, x[1], x[2] * 1e5))
    g = ForwardDiff.gradient(lnrho, [T, P_bar])
    H = ForwardDiff.hessian(lnrho, [T, P_bar])
    alpha = -g[1]
    beta = g[2]
    daldT = -H[1, 1]
    rho = MW_H2O / V0
    return water_from_density(T, rho, alpha, daldT, beta)
end

# ---------------------------------------------------------------------------
# g-function (Shock et al. 1992, Eqs. 24-26, 32-33), as in SUPCRT92/CHNOSZ.
# Only g itself is needed for G. D in g/cm3, Tc in C, P in bar.
# ---------------------------------------------------------------------------
function gfun(D, Tc, P)
    D >= 1 && return 0.0
    P > 6000 && return 0.0
    ag = -2.037662 + 5.747000e-3 * Tc - 6.557892e-6 * Tc^2
    bg = 6.107361 - 1.074377e-2 * Tc + 1.268348e-5 * Tc^2
    g = ag * (1 - D)^bg
    if Tc > 155 && P < 1000 && Tc < 355
        u = (Tc - 155) / 300
        f = (u^4.8 + 0.3666666e2 * u^16) * (-0.1504956e-9 * (1000 - P)^3 + 0.5017997e-13 * (1000 - P)^4)
        g -= f
    end
    return g
end

struct HKFRecord
    name::String
    Z::Float64
    Gf::Float64   # cal/mol (apparent Gf at Tr, Pr)
    Hf::Float64
    Sr::Float64   # cal/mol/K
    a1::Float64   # cal/mol/bar
    a2::Float64   # cal/mol
    a3::Float64   # cal K/mol/bar
    a4::Float64   # cal K/mol
    c1::Float64   # cal/mol/K
    c2::Float64   # cal K/mol
    w::Float64    # omega(Tr, Pr), cal/mol
end

# supcrt92.xml stores the HKF parameters in their true units (as read by the app);
# standard_states.hkf_ion_g0 scales them to CHNOSZ units and CHNOSZ scales back.
hkf_records_from_dict(d) = HKFRecord(String(d["Species"]), d["Charge"], d["Gf"], d["Hf"], d["Sr"],
                                     d["a1"], d["a2"], d["a3"], d["a4"], d["c1"], d["c2"], d["wref"])

"""Apparent standard Gibbs energy (cal/mol, CHNOSZ convention) of an aqueous
species at T [K], P [bar] given the water state there and at Tr, Pr."""
function hkf_G(r::HKFRecord, T, P, w::WaterState, wr::WaterState)
    # nonsolvation
    lnTr = log((TR * (T - THETA)) / (T * (TR - THETA)))
    Gc = -r.c1 * (T * log(T / TR) - T + TR) -
         r.c2 * ((1 / (T - THETA) - 1 / (TR - THETA)) * ((THETA - T) / THETA) - (T / THETA^2) * lnTr)
    lnP = log((PSI + P) / (PSI + PR))
    Ga = r.a1 * (P - PR) + r.a2 * lnP + (r.a3 * (P - PR) + r.a4 * lnP) / (T - THETA)
    # solvation with omega(T, P) from the g-function (SUPCRT92 reac92.f)
    Z = r.Z
    omega = r.w
    omegaPT = omega
    if Z != 0
        g = gfun(w.rho / 1000, T - 273.15, P)
        reref = Z^2 / (omega / ETA + Z / 3.082)
        re = reref + abs(Z) * g
        omegaPT = ETA * (Z^2 / re - Z / (3.082 + g))
    end
    ZB = -1 / w.epsilon
    ZBr = -1 / wr.epsilon
    Gs = -omegaPT * (ZB + 1) + omega * (ZBr + 1) + omega * wr.YBorn * (T - TR)
    # origination
    Go = r.Gf - r.Sr * (T - TR)
    return Gc + Ga + Gs + Go
end

"""The app's mu0 (J/mol, NASA Hf - T S reference) from the apparent HKF G."""
app_mu0(r::HKFRecord, G_cal) = (G_cal - (r.Gf - r.Hf + TR * r.Sr)) * 4.184

end # module
