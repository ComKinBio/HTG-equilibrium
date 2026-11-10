"""Standard-state Gibbs energies g°(T, P) [J/mol].

All sources share one reference: Hf(298.15 K) - T S, as used by NASA-7.

    molecules       NASA-7 polynomials (pressure independent, ideal-gas standard state)
    NASA solids     NASA-7 polynomials
    SUPCRT minerals Maier-Kelley Cp = a + bT + c/T^2 plus a constant-volume term Vr (P - Pr)
    aqueous ions    HKF (CHNOSZ), shifted from CHNOSZ's apparent Gf(298.15 K) anchor
                    to the NASA reference
"""
import numpy as np
import pandas as pd

from .backends import chnosz
from .constants import CAL_TO_J, P_REF, T_REF
from .databases import nasa7_properties


def nasa_g0(registry, nasa_db, species, T):
    """g° of fluid molecules (or NASA solids) from their NASA-7 records."""
    mu0 = np.zeros(len(species), dtype=float)
    for i, name in enumerate(species):
        nasa_name = registry.nasa_name(name)
        if nasa_name not in nasa_db:
            raise KeyError(f"NASA species '{nasa_name}' for '{name}' is not present in nasa_db.")
        mu0[i] = float(nasa7_properties(nasa_db, nasa_name, T)[3])
    return mu0


def supcrt_mineral_properties(params, name, T, P):
    """Mineral G, H, S, Cp, V from the single-range Maier-Kelley record ``params``.

    H(T) = Hf(Tr) + int Cp dT,  S(T) = Sr + int Cp/T dT,  G = H - T S + Vr (P - Pr).
    Tmax is the equation limit, not a hydrate stability temperature.
    """
    t, tr, p = float(T), T_REF, float(P)
    if not np.isfinite(t) or not tr <= t <= params['Tmax']:
        raise ValueError(f"{name}: T={t:g} K outside supported [{tr:g}, {params['Tmax']:g}] K")
    if not np.isfinite(p) or p <= 0:
        raise ValueError('Invalid mineral pressure')
    a, b, c = [params[k] * CAL_TO_J for k in ('a', 'b', 'c')]
    dh = a * (t - tr) + b * (t * t - tr * tr) / 2 - c * (1 / t - 1 / tr)
    ds = a * np.log(t / tr) + b * (t - tr) - c * (t ** -2 - tr ** -2) / 2
    h, entropy = params['Hf'] * CAL_TO_J + dh, params['Sr'] * CAL_TO_J + ds
    pv = params['Vr'] * 1e-6 * (p - P_REF)
    return {'G': h - t * entropy + pv, 'H': h + pv, 'S': entropy,
            'Cp': a + b * t + c / t ** 2, 'V_m3_per_mol': params['Vr'] * 1e-6}


def solid_g0(registry, nasa_db, species, T, P):
    """g° of pure solids: SUPCRT minerals or NASA-7 solids."""
    values = []
    for name in species:
        d = registry.solid[name]
        if d.get('adapter') == 'supcrt_mineral':
            value = supcrt_mineral_properties(d['supcrt_mineral'], name, T, P)['G']
        else:
            value = nasa7_properties(nasa_db, registry.nasa_name(name), T)[3]
        values.append(value)
    return np.asarray(values)


def supcrt_water(T, P):
    """SUPCRT92 water properties at T [K], P [Pa] (density in kg/m3, permittivity, ...)."""
    props = ['rho', 'epsilon', 'QBorn', 'XBorn', 'YBorn', 'alpha', 'daldT', 'beta']
    out = chnosz().water(props, T=float(T), P=float(P) / 100000.0)
    return {key: float(np.asarray(out[key]).reshape(-1)[0]) for key in props}


def hkf_ion_g0(hkf_table, ions, T, P):
    """g° of aqueous ions from their HKF records, on the NASA (Hf - T S) reference."""
    records = hkf_table.set_index('Species')
    rows = []
    for species in ions:
        d = records.loc[species]
        rows.append({'name': str(d['SUPCRT name']), 'state': 'aq', 'G': float(d['Gf']),
                     'H': float(d['Hf']), 'S': float(d['Sr']), 'Cp': np.nan, 'V': np.nan,
                     'a1.a': float(d['a1']) * 10, 'a2.b': float(d['a2']) / 100,
                     'a3.c': float(d['a3']), 'a4.d': float(d['a4']) / 10000.0,
                     'c1.e': float(d['c1']), 'c2.f': float(d['c2']) / 10000.0,
                     'omega.lambda': float(d['wref']) / 100000.0, 'z.T': float(d['Charge'])})
    params = pd.DataFrame(rows)
    result, _ = chnosz().hkf(property=['G'], parameters=params, T=float(T), P=float(P) / 100000.0,
                             contrib=['n', 's', 'o'], H2O_props=['rho'], water_model='SUPCRT92')
    values = np.array([float(np.asarray(result[i]['G']).reshape(-1)[0]) * 4.184
                       for i in range(len(rows))])
    if not np.all(np.isfinite(values)):
        raise ValueError('Nonfinite HKF ion standard Gibbs energies')
    # CHNOSZ returns the apparent G anchored to Gf(Tr); NASA uses Hf(Tr) - Tr*Sr.
    # Shift the constant reference of each ion. Never subtract T*Sr here: HKF
    # already includes the entropy evolution.
    offsets = np.array([(float(records.loc[s, 'Gf']) - float(records.loc[s, 'Hf'])
                         + 298.15 * float(records.loc[s, 'Sr'])) * 4.184 for s in ions])
    return values - offsets
