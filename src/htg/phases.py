"""Phase models: Gibbs energy G(n) and chemical potentials mu(n) of one phase.

Every phase is built for fixed T [K] and P [Pa] and exposes

    gibbs(n)          G [J] of amounts n [mol], ordered like ``species``
    gibbs_and_mu(n)   (G, mu): mu_i = dG/dn_i; G = sum n_i mu_i (Euler)
    hessian(n)        dmu/dn (fluid phases only)

Formulas
--------
Fluid molecules (PC-SAFT):  mu_i = g°_i(T) + RT ln(x_i P phi_i / P0)
Pure solids:                mu_i = g°_i(T, P)
Aqueous electrolyte (ePC-SAFT + HKF):
    neutrals  mu_i = g°_i(T) + RT ln(x_i P phi_i / P0)
    ions      mu_i = g°_i,HKF(T, P) + RT [ln m_i + ln x_w + ln phi_i - ln phi_i^inf]
with molality m_i = n_i / (M_w n_w) and phi_i^inf the infinite-dilution limit.
"""
import numpy as np
import pandas as pd

from . import standard_states as ss
from .backends import clapeyron, julia, julia_int_vector, julia_string_vector
from .constants import ATOMIC_WEIGHTS, MW_WATER_KG, P0, R
from .species import molecular_weight


class EOSFailure(ValueError):
    """The equation of state gave no valid fugacity coefficients at this composition."""


def association_sites(association_scheme, epsilonAB, kappaAB):
    scheme = str(association_scheme).strip().upper()
    if scheme == '2B':
        return 1, 1
    if scheme == '3B':
        return 1, 2
    if np.isfinite(epsilonAB) and np.isfinite(kappaAB) and epsilonAB > 0.0 and kappaAB > 0.0:
        return 1, 1
    return 0, 0


def _first_value(mapping, names, default=np.nan):
    for name in names:
        if name in mapping:
            value = mapping[name]
            if value is not None and not pd.isna(value):
                return value
    return default


def _required_float(value, parameter, species):
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{species}: invalid {parameter} = {value!r}') from exc
    if not np.isfinite(result):
        raise ValueError(f'{species}: missing {parameter}')
    return result


def _as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def build_pcsaft_model(species, resolved, binary_interactions):
    """Clapeyron PCSAFT model for molecules ``species`` with resolved parameter sets."""
    components, MW, segment, sigma, epsilon, n_H, n_e = [], [], [], [], [], [], []
    epsilon_assoc, bondvol = {}, {}
    for name in species:
        if name not in resolved:
            raise KeyError(f"No resolved PC-SAFT parameter set for '{name}'.")
        p = resolved[name]
        component = str(p['Formula']).strip()
        components.append(component)
        MW.append(molecular_weight(component, ATOMIC_WEIGHTS))
        segment.append(float(p['m']))
        sigma.append(float(p['sigma']))
        epsilon.append(float(p['epsilon']))
        epsilonAB = _as_float(p.get('epsilonAB', np.nan))
        kappaAB = _as_float(p.get('kappaAB', np.nan))
        nH, ne = association_sites(p.get('association_scheme', '-'), epsilonAB, kappaAB)
        n_H.append(nH)
        n_e.append(ne)
        if nH > 0 and ne > 0 and np.isfinite(epsilonAB) and np.isfinite(kappaAB):
            key = ((component, 'e'), (component, 'H'))
            epsilon_assoc[key] = float(epsilonAB)
            bondvol[key] = float(kappaAB)
    MW, segment = np.asarray(MW, dtype=float), np.asarray(segment, dtype=float)
    sigma, epsilon = np.asarray(sigma, dtype=float), np.asarray(epsilon, dtype=float)
    n_H, n_e = np.asarray(n_H, dtype=int), np.asarray(n_e, dtype=int)
    for label, values in {'Mw': MW, 'segment': segment, 'sigma': sigma, 'epsilon': epsilon}.items():
        if np.any(~np.isfinite(values)):
            raise ValueError(f'{label} contains non-finite values.')
    userlocations = {'Mw': MW, 'segment': segment, 'sigma': sigma, 'epsilon': epsilon}
    if np.any(n_H > 0) or np.any(n_e > 0):
        userlocations['n_H'] = n_H
        userlocations['n_e'] = n_e
        if epsilon_assoc:
            userlocations['epsilon_assoc'] = epsilon_assoc
        if bondvol:
            userlocations['bondvol'] = bondvol
    k_matrix = None
    for (a, b), kij in binary_interactions:
        if a in components and b in components:
            if k_matrix is None:
                k_matrix = np.zeros((len(components), len(components)))
            i, j = components.index(a), components.index(b)
            k_matrix[i, j] = k_matrix[j, i] = kij
    if k_matrix is not None:
        userlocations['k'] = k_matrix
    return clapeyron().PCSAFT(components, userlocations=userlocations)


class SolidPhase:
    """Pure solids: G = sum n_i g°_i(T, P)."""

    def __init__(self, species, mu0):
        self.species = list(species)
        self.mu0 = np.asarray(mu0, dtype=float)

    def gibbs(self, n):
        return float(n @ self.mu0)

    def gibbs_and_mu(self, n):
        return float(n @ self.mu0), self.mu0


class FluidPhase:
    """Molecular fluid phase with PC-SAFT fugacity coefficients.

    root: 'vapour' (gas) or 'liquid' (organic liquid, molecular aqueous). The root
    choice alone does not establish phase stability.
    """

    def __init__(self, name, species, mu0, eos_model, T, P, root):
        self.name = name
        self.species = list(species)
        self.mu0 = np.asarray(mu0, dtype=float)
        self.eos_model = eos_model
        self.T, self.P, self.root = float(T), float(P), root
        self._pure_mu = None

    def _phis(self, positive):
        phis = np.asarray(clapeyron().fugacity_coefficient(self.eos_model, self.P, self.T, positive,
                                                           phase=self.root),
                          dtype=float).reshape(-1)
        if np.any(~np.isfinite(phis)) or np.any(phis <= 0):
            raise EOSFailure(f'Invalid PC-SAFT fugacities: phase={self.root}, '
                             f'T={self.T - 273.15:g} C, P={self.P / 100000.0:g} bar, '
                             f'species={self.species}, x={(positive / positive.sum()).tolist()}')
        return phis

    def mu_at(self, n):
        positive = np.maximum(n, 1e-14)
        phis = self._phis(positive)
        return self.mu0 + R * self.T * np.log(positive / positive.sum() * self.P / P0 * phis)

    def gibbs(self, n):
        if np.sum(n) <= 1e-13:
            return 0.0
        return float(np.dot(n, self.mu_at(n)))

    def gibbs_and_mu(self, n):
        if np.sum(n) <= 1e-13:
            # G = 0 for an empty phase. G is homogeneous of degree one, so the
            # one-sided derivative along e_i is the pure-i chemical potential.
            if self._pure_mu is None:
                eye = np.eye(len(self.species))
                self._pure_mu = np.array([self.mu_at(eye[i])[i] for i in range(len(self.species))])
            return 0.0, self._pure_mu.copy()
        mu = self.mu_at(n)
        return float(n @ mu), mu

    def lnphi(self, n):
        return np.log(self._phis(np.maximum(n, 1e-14)))

    def hessian(self, n):
        if np.sum(n) <= 1e-13:
            return np.zeros((len(n), len(n)))
        return _hessian(self, n, np.maximum(n, 1e-14).sum())


class ElectrolytePhase:
    """Aqueous electrolyte: ePC-SAFT (PC-SAFT + Debye-Hückel) with HKF ion standard states.

    The equation of state always contains the registry's reference neutrals and
    ions (in that order) followed by the selected species; unselected
    components carry trace probe amounts. ``species`` lists the case's selected
    aqueous species; ``all_species`` the full component vector.
    """

    def __init__(self, model, selected, T, P, resolved, atomic_weights):
        registry = model.registry
        ions = registry.ion
        self.T, self.P = float(T), float(P)
        self.species = list(selected)
        self.neutrals = list(dict.fromkeys(registry.aqueous_reference_neutrals
                                           + [s for s in selected if s not in ions]))
        self.ions = list(dict.fromkeys(registry.aqueous_reference_ions
                                       + [s for s in selected if s in ions]))
        self.all_species = self.neutrals + self.ions
        if self.all_species[0] != 'water':
            raise ValueError('The aqueous reference species must start with water.')
        self.positions = np.array([self.all_species.index(s) for s in selected], dtype=int)
        self.z = np.array([ions[s]['charge'] for s in self.ions], dtype=float)
        self.probe_weights = registry.probe_weights

        self.water = ss.supcrt_water(self.T, self.P)
        eps = self.water['epsilon']
        if not (np.isfinite(self.water['rho']) and self.water['rho'] >= 350
                and np.isfinite(eps) and eps > 1):
            raise ValueError(f"State outside SUPCRT aqueous validity checks: "
                             f"rho={self.water['rho']}, epsilon={eps}")
        self.short_range = self._build_short_range(model, resolved, atomic_weights)
        cl = clapeyron()
        refs = julia().seval('String[]')
        components = julia_string_vector(self.components)
        ion_model = cl.hsdDH(components, cl.ConstRSP(float(eps)), refs)
        self.eos_model = cl.ePCSAFT(components, julia_int_vector(self.charges),
                                    self.short_range.idealmodel, self.short_range, ion_model, refs)
        self.mu0_ion = ss.hkf_ion_g0(registry.ion_hkf_table, self.ions, self.T, self.P)
        self.mu0_neutral = ss.nasa_g0(registry, model.databases.nasa, self.neutrals, self.T)
        self.lnphi_infinity = self._infinite_dilution_reference()

    def _build_short_range(self, model, resolved, atomic_weights):
        registry = model.registry
        all_molecular = registry.all_molecular
        self.components = ([all_molecular[s]['formula'] for s in self.neutrals]
                           + [registry.ion[s]['formula'] for s in self.ions])
        self.charges = [0] * len(self.neutrals) + [int(registry.ion_charge[s]) for s in self.ions]
        ion_lookup = {str(row['Species']): row.to_dict()
                      for _, row in registry.ion_epcsaft_table.iterrows()}
        ion_mw = {name: sum(atomic_weights[e] * count for e, count in d['elements'].items())
                  for name, d in registry.ion.items()}
        Mw, segment, sigma, epsilon, n_H, n_e = [], [], [], [], [], []
        epsilon_assoc, bondvol = {}, {}
        for index, name in enumerate(self.all_species):
            if name in self.neutrals:
                if name not in resolved:
                    raise KeyError(f'No resolved molecular PC-SAFT parameters for {name}')
                p = resolved[name]
                default_mw = sum(atomic_weights[e] * c for e, c in all_molecular[name]['elements'].items())
            else:
                if name not in ion_lookup:
                    raise KeyError(f'No resolved ion ePC-SAFT parameters for {name}')
                p = ion_lookup[name]
                default_mw = ion_mw[name]
            Mw_i = _required_float(_first_value(p, ['Mw', 'MW', 'MW_g_mol', 'MW [g/mol]'],
                                                default=default_mw), 'Mw', name)
            m_i = _required_float(_first_value(p, ['m', 'segment']), 'm', name)
            sigma_i = _required_float(_first_value(p, ['sigma', 'sigma_A', 'sigma [A]']), 'sigma', name)
            eps_i = _required_float(_first_value(p, ['epsilon', 'u_over_kB_K', 'epsilon/kB [K]']),
                                    'epsilon', name)
            epsAB = _first_value(p, ['epsilonAB', 'epsilon_assoc', 'epsilonAB/kB [K]'])
            kapAB = _first_value(p, ['kappaAB', 'bondvol'])
            scheme = _first_value(p, ['association_scheme', 'Association'], default='-')
            Mw.append(Mw_i)
            segment.append(m_i)
            sigma.append(sigma_i)
            epsilon.append(eps_i)
            if (pd.notna(epsAB) and pd.notna(kapAB) and np.isfinite(float(epsAB))
                    and np.isfinite(float(kapAB))):
                epsAB, kapAB = float(epsAB), float(kapAB)
                nH, ne = association_sites(scheme, epsAB, kapAB)
            else:
                nH, ne = 0, 0
            n_H.append(int(nH))
            n_e.append(int(ne))
            component = self.components[index]
            if nH > 0 and ne > 0:
                key = ((component, 'H'), (component, 'e'))
                epsilon_assoc[key] = epsAB
                bondvol[key] = float(kapAB)
        userlocations = {'Mw': np.asarray(Mw, dtype=float), 'segment': np.asarray(segment, dtype=float),
                         'sigma': np.asarray(sigma, dtype=float),
                         'epsilon': np.asarray(epsilon, dtype=float),
                         'n_H': np.asarray(n_H, dtype=int), 'n_e': np.asarray(n_e, dtype=int),
                         'epsilon_assoc': epsilon_assoc, 'bondvol': bondvol}
        return clapeyron().PCSAFT(self.components, userlocations=userlocations, verbose=False)

    # ---------------------------------------------------------------- EOS probes
    def _balanced_probe_weights(self):
        weights = np.array([self.probe_weights.get(s, 1.) for s in self.ions], dtype=float)
        imbalance = float(self.z @ weights)
        candidates = np.flatnonzero(self.z * imbalance < 0)
        if abs(imbalance) > 1e-14:
            if not candidates.size:
                raise ValueError('Both positive and negative ionic carriers are needed for EOS probes')
            j = candidates[0]
            weights[j] += -imbalance / self.z[j]
        return weights

    def phi(self, n):
        """Fugacity coefficients for the full component vector ``n``.

        The probe composition is made electroneutral (smoothly) and every component
        gets a trace amount; the physical amounts used in balances are unchanged.
        """
        n = np.asarray(n, dtype=float).copy()
        if n.shape != (len(self.all_species),) or np.any(n < 0):
            raise ValueError('Invalid aqueous composition vector')
        split, z = len(self.neutrals), self.z
        imbalance = float(n[split:] @ z)
        smoothing = np.hypot(imbalance, 2e-14)
        pos, neg = np.flatnonzero(z > 0), np.flatnonzero(z < 0)
        if not pos.size or not neg.size:
            raise ValueError('Aqueous EOS probe requires positive and negative ionic carriers')
        n[split + pos[0]] += (smoothing - imbalance) / (2 * z[pos[0]])
        n[split + neg[0]] += (smoothing + imbalance) / (-2 * z[neg[0]])
        n[:split] += 1e-14
        n[split:] += self._balanced_probe_weights() * 1e-14
        phi = np.asarray(clapeyron().fugacity_coefficient(self.eos_model, self.P, self.T, n,
                                                          phase='liquid'), dtype=float).reshape(-1)
        if len(phi) != len(n) or np.any(~np.isfinite(phi)) or np.any(phi <= 0):
            raise EOSFailure('Invalid ePC-SAFT fugacity coefficients')
        return phi

    def _infinite_dilution_reference(self):
        logs = []
        split = len(self.neutrals)
        for scale in (1e-5, 1e-6, 1e-7):
            n = np.zeros(len(self.all_species))
            n[self.all_species.index('water')] = 10.
            n[split:] = scale * self._balanced_probe_weights()
            logs.append(np.log(self.phi(n)))
        x = np.sqrt(np.array([1e-5, 1e-6, 1e-7]))
        intercept = np.array([np.polyfit(x, np.array(logs)[:, j], 1)[1]
                              for j in range(len(self.all_species))])
        estimate = logs[-1][split:] - (logs[-2][split:] - logs[-1][split:]) * x[-1] / (x[-2] - x[-1])
        if np.max(np.abs(intercept[split:] - estimate)) > .02:
            raise ValueError('Ion dilution extrapolation did not converge')
        return intercept

    # ---------------------------------------------------------------- G and mu
    def expand(self, n):
        full = np.zeros(len(self.all_species))
        full[self.positions] = n
        return full

    def gibbs_full(self, n):
        """G of the full component vector (water first)."""
        if n[0] <= 1e-09:
            if np.sum(n[1:]) > 1e-10:
                raise ValueError('Aqueous solutes without water')
            return 0.0
        phi = self.phi(n)
        q = np.maximum(n, 1e-14)
        x = q / q.sum()
        xw = x[0]
        k, s = len(self.ions), len(self.neutrals)
        mu_neutral = self.mu0_neutral + R * self.T * np.log(x[:s] * self.P / P0 * phi[:s])
        molal = np.maximum(n[-k:], 1e-14) / (MW_WATER_KG * n[0])
        mu_ions = self.mu0_ion + R * self.T * (np.log(molal) + np.log(xw) + np.log(phi[-k:])
                                               - self.lnphi_infinity[-k:])
        return float(np.dot(n[:s], mu_neutral) + np.dot(n[-k:], mu_ions))

    def gibbs(self, n):
        return self.gibbs_full(self.expand(n))

    def _phi_for_optimizer(self, full):
        # On the optimizer path every probe failure is an EOS failure (cut the step back).
        try:
            return self.phi(full)
        except ValueError as exc:
            raise EOSFailure(str(exc)) from exc

    def gibbs_and_mu(self, n):
        full = self.expand(n)
        if full[0] <= 1e-9:
            if np.sum(full[1:]) > 1e-10:
                raise ValueError('Aqueous solutes without water')
            return 0.0, np.zeros(len(n))
        phi = self._phi_for_optimizer(full)
        RT, split, n_ion = R * self.T, len(self.neutrals), len(self.ions)
        q = np.maximum(full, 1e-14)
        x = q / q.sum()
        mu_neutral = self.mu0_neutral + RT * np.log(x[:split] * self.P / P0 * phi[:split])
        molal = np.maximum(full[-n_ion:], 1e-14) / (MW_WATER_KG * full[0])
        mu_ions = self.mu0_ion + RT * (np.log(molal) + np.log(x[0]) + np.log(phi[-n_ion:])
                                       - self.lnphi_infinity[-n_ion:])
        G = float(full[:split] @ mu_neutral + full[-n_ion:] @ mu_ions)
        return G, np.concatenate([mu_neutral, mu_ions])[self.positions]

    def mu_at(self, n):
        return self.gibbs_and_mu(n)[1]

    def lnphi(self, n):
        return np.log(self._phi_for_optimizer(self.expand(n)))[self.positions]

    def hessian(self, n):
        if n[0] <= 1e-9:
            return np.zeros((len(n), len(n)))
        # Unselected reference components enter the mole fractions at 1e-14.
        total = np.maximum(n, 1e-14).sum() + 1e-14 * (len(self.all_species) - len(n))
        return _hessian(self, n, total)


def _hessian(phase, n, total):
    """d mu / d n: analytic ideal part + forward-difference ln(phi) part, symmetrized."""
    RT = R * phase.T
    alive = n > 1e-14
    ideal = np.diag(np.where(alive, 1.0 / np.maximum(n, 1e-14), 0.0)) - np.outer(alive, alive) / total
    base = phase.lnphi(n)
    jac = np.zeros((len(n), len(n)))
    h = 1e-7 * n.sum()
    for j in range(len(n)):
        step = n.copy()
        step[j] += h
        jac[:, j] = (phase.lnphi(step) - base) / h
    return RT * (ideal + 0.5 * (jac + jac.T))
