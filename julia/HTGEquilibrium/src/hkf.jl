# Revised HKF standard Gibbs energy of aqueous ions (Tanger & Helgeson 1988;
# Shock et al. 1992), with omega(T, P) from the g-function. Python
# counterpart: `standard_states.hkf_ion_g0`, which calls pychnosz `hkf`.
# The arithmetic follows pychnosz (models/hkf.py) term by term, including the
# parameter scaling round trip of `hkf_ion_g0`, so that the results agree to
# ~1e-13 rather than only to the reference tolerance.

const THETA = 228.0        # K
const PSI = 2600.0         # bar
const ETA = 1.66027e5      # Angstrom cal/mol

"""
HKF record of one aqueous species, in the units of `data/supcrt92.xml` as the
app's registry reads them (`registry.ion_hkf_table`). The constructor from a
dict applies the same scaling round trip as Python `hkf_ion_g0` and pychnosz
(e.g. a1 -> a1*10 -> *10^-1), so the Float64 values are identical.
"""
struct HKFRecord
    species::String       # app name, e.g. "CO3^2-"
    supcrt_name::String   # SUPCRT92 name, e.g. "CO3-2"
    Z::Float64            # charge
    Gf::Float64           # cal/mol, apparent Gf at Tr, Pr
    Hf::Float64           # cal/mol
    Sr::Float64           # cal/(mol K)
    a1::Float64           # cal/(mol bar)
    a2::Float64           # cal/mol
    a3::Float64           # cal K/(mol bar)
    a4::Float64           # cal K/mol
    c1::Float64           # cal/(mol K)
    c2::Float64           # cal K/mol
    omega::Float64        # cal/mol, omega at Tr, Pr
end

"""
    HKFRecord(d)

From a dict with the columns of the app's `ion_hkf_table`: `Species`,
`SUPCRT name`, `Charge`, `Gf`, `Hf`, `Sr`, `a1`, `a2`, `a3`, `a4`, `c1`, `c2`,
`wref`.
"""
function HKFRecord(d::AbstractDict)
    f(k) = Float64(d[k])
    return HKFRecord(String(d["Species"]), String(d["SUPCRT name"]), f("Charge"), f("Gf"), f("Hf"), f("Sr"),
                     (f("a1") * 10) * 10.0^-1, (f("a2") / 100) * 10.0^2, f("a3"),
                     (f("a4") / 10000.0) * 10.0^4, f("c1"), (f("c2") / 10000.0) * 10.0^4,
                     (f("wref") / 100000.0) * 10.0^5)
end

"""
    gfun(D, Tc, P)

g-function of Shock et al. (1992), Eqs. 24-26 and 32-33, as in SUPCRT92 and
CHNOSZ. D: water density [g/cm3], Tc: temperature [°C], P: pressure [bar].
"""
function gfun(D, Tc, P)
    D < 1 || return 0.0
    ag = -2.037662 + 5.747000e-3 * Tc + -6.557892e-6 * Tc^2
    bg = 6.107361 + -1.074377e-2 * Tc + 1.268348e-5 * Tc^2
    g = ag * (1 - D)^bg
    if Tc > 155 && P < 1000 && Tc < 355
        u = (Tc - 155) / 300
        f = (u^4.8 + 0.3666666e2 * u^16) * (-0.1504956e-9 * (1000 - P)^3 + 0.5017997e-13 * (1000 - P)^4)
        g -= f
    end
    P > 6000 && return 0.0
    return g
end

"""
    hkf_apparent_G(r::HKFRecord, T, P_bar, w, wr) -> cal/mol

Apparent standard Gibbs energy (CHNOSZ convention, anchored at Gf(Tr, Pr)) at
T [K] and P [bar], with the water properties `w` at (T, P) and `wr` at
(25 °C, 1 bar). Sum of the nonsolvation, solvation and origination terms.
"""
function hkf_apparent_G(r::HKFRecord, T, P, w::WaterProperties, wr::WaterProperties)
    at_ref = T == TR && P == PR
    # omega(T, P) after SUPCRT92 reac92.f; H+ keeps omega = 0 by convention.
    omega = r.omega
    omega_PT = omega
    Z = r.Z
    if Z != 0 && r.supcrt_name != "H+"
        g = gfun(w.rho / 1000, T - 273.15, P)
        reref = Z^2 / (omega / ETA + Z / (3.082 + 0))
        re = reref + abs(Z) * g
        omega_PT = ETA * (Z^2 / re - Z / (3.082 + g))
    end
    # nonsolvation
    lnP = log((PSI + P) / (PSI + PR))
    p_c = -r.c1 * (T * log(T / TR) - T + TR) -
          r.c2 * ((1 / (T - THETA) - 1 / (TR - THETA)) * ((THETA - T) / THETA) -
                  (T / THETA^2) * log((TR * (T - THETA)) / (T * (TR - THETA))))
    p_a = r.a1 * (P - PR) + r.a2 * lnP + (r.a3 * (P - PR) + r.a4 * lnP) / (T - THETA)
    n = at_ref ? 0.0 : p_c + p_a
    # solvation
    ZB = -1 / w.epsilon
    ZBr = -1 / wr.epsilon
    s = -omega_PT * (ZB + 1) + omega * (ZBr + 1) + omega * wr.YBorn * (T - TR)
    s = at_ref ? 0.0 : s
    # origination
    o = T == TR ? r.Gf : r.Gf - r.Sr * (T - TR)
    return ((0.0 + n) + s) + o
end

"""
    hkf_ion_g0(model::WaterModel, records, T, P) -> Vector{Float64}

Standard chemical potentials [J/mol] of aqueous species at T [K], P [Pa] on
the app's reference (NASA Hf(298.15 K) - T S), as Python `hkf_ion_g0`: the
apparent HKF G (cal/mol) times 4.184, minus the constant shift
(Gf - Hf + 298.15 Sr) * 4.184. All water properties, including the 25 °C,
1 bar reference state, come from `model`: one water model is used throughout,
never a mix.
"""
function hkf_ion_g0(model::WaterModel, records::AbstractVector{HKFRecord}, T::Real, P::Real)
    reference = water_properties(model, TR, PR_PA)
    T = Float64(T)
    P_bar = Float64(P) / 100000.0
    w = water_properties(model, T, P)
    values = [hkf_apparent_G(r, T, P_bar, w, reference) * 4.184 for r in records]
    all(isfinite, values) || throw(DomainError((T, P), "Nonfinite HKF ion standard Gibbs energies"))
    offsets = [(r.Gf - r.Hf + 298.15 * r.Sr) * 4.184 for r in records]
    return values .- offsets
end
