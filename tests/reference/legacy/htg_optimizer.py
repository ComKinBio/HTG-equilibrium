"""Optimization layer for HTG_Model.solve_equilibrium.

The Gibbs model (EOS, standard states, phase definitions) belongs to the model
notebook; this module only decides how G is minimized.  It evaluates G and its
exact gradient.  For every phase G = sum_i n_i mu_i with mu_i from the same
Helmholtz model, so dG/dn_i = mu_i (Euler/Gibbs-Duhem): no AD or finite
differences are needed.  The chemical potentials repeat the formulas of
``molecular_gibbs`` and ``aqueous_gibbs``; ``GibbsProblem.check_against_model``
asserts that both give the same G, so a change to the model is detected rather
than silently diverging.
"""
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, minimize

EOS_FAILURE_TEXT = ('Invalid PC-SAFT fugacities', 'Invalid ePC-SAFT fugacity coefficients',
                    'Nonfinite phase Gibbs energy')


class EOSFailure(ValueError):
    pass


class GibbsProblem:
    """min_x (G(expand(x)) - G0)/scale  s.t.  A x = b,  c x >= d,  x >= lower.

    ``x`` holds the free (phase, species) amounts; ``ns`` is the model module's
    global namespace, read at call time because solve_equilibrium swaps T, P and
    the aqueous state in and out of it.
    """

    def __init__(self, ns, pairs, free, phase_species, phase_indices, mu0, aqueous_model,
                 aq_positions, A_eq, b_eq, ineq_row, ineq_target, lower, x_start, G0, scale,
                 model_total_G):
        self.ns = ns
        self.pairs = list(pairs)
        self.free = np.asarray(free, dtype=int)
        self.phase_species = phase_species
        self.phase_indices = phase_indices
        self.mu0 = mu0
        self.aqueous_model = aqueous_model
        self.aq_positions = aq_positions
        self.A_eq = np.asarray(A_eq, dtype=float)
        self.b_eq = np.asarray(b_eq, dtype=float)
        self.ineq_row = None if ineq_row is None else np.asarray(ineq_row, dtype=float)
        self.ineq_target = ineq_target
        self.lower = np.asarray(lower, dtype=float)
        self.x_start = np.asarray(x_start, dtype=float)
        self.G0 = float(G0)
        self.scale = float(scale)
        self.model_total_G = model_total_G
        self.eos_errors = []
        self.nfev = 0
        self._eos_models = {}
        self._pure_mu = {}
        self._last = None

    # ----- chemical potentials (same formulas as the model notebook) -----
    def _eos_model(self, names):
        key = tuple(names)
        if key not in self._eos_models:
            self._eos_models[key], _ = self.ns['build_pcsaft_model'](list(names))
        return self._eos_models[key]

    def _molecular_mu_at(self, n, names, mu0, root):
        ns = self.ns
        positive = np.maximum(n, 1e-14)
        phis = np.asarray(ns['cl'].fugacity_coefficient(self._eos_model(names), ns['P'], ns['T'],
                                                        positive, phase=root), dtype=float).reshape(-1)
        if np.any(~np.isfinite(phis)) or np.any(phis <= 0):
            raise EOSFailure(f'Invalid PC-SAFT fugacities: phase={root}, species={list(names)}, '
                             f'x={(positive / positive.sum()).tolist()}')
        return mu0 + ns['R'] * ns['T'] * np.log(positive / positive.sum() * ns['P'] / ns['P0'] * phis)

    def _molecular(self, phase, n, root):
        names, mu0 = self.phase_species[phase], self.mu0[phase]
        if np.sum(n) <= 1e-13:
            # The model sets G = 0 for an empty phase.  G is homogeneous of degree
            # one, so the one-sided derivative along e_i is the pure-i potential.
            if phase not in self._pure_mu:
                self._pure_mu[phase] = np.array([
                    self._molecular_mu_at(np.eye(len(names))[i], names, mu0, root)[i]
                    for i in range(len(names))])
            return 0.0, self._pure_mu[phase].copy()
        mu = self._molecular_mu_at(n, names, mu0, root)
        return float(n @ mu), mu

    def _electrolyte(self, n):
        ns = self.ns
        full = np.zeros(len(ns['AQUEOUS_SPECIES']))
        full[self.aq_positions] = n
        if full[0] <= 1e-9:
            if np.sum(full[1:]) > 1e-10:
                raise ValueError('Aqueous solutes without water')
            return 0.0, np.zeros(len(n))
        try:
            phi = ns['phi_electrolyte'](full)
        except ValueError as exc:
            raise EOSFailure(str(exc)) from exc
        RT, split = ns['R'] * ns['T'], len(ns['AQ_NEUTRAL_SPECIES'])
        n_ion = len(ns['AQ_ION_SPECIES'])
        q = np.maximum(full, 1e-14)
        x = q / q.sum()
        mu_neutral = ns['mu0_aq_neutral'] + RT * np.log(x[:split] * ns['P'] / ns['P0'] * phi[:split])
        molal = np.maximum(full[-n_ion:], 1e-14) / (ns['MW_WATER_KG'] * full[0])
        mu_ions = ns['mu0_ion'] + RT * (np.log(molal) + np.log(x[0]) + np.log(phi[-n_ion:])
                                        - ns['lnphi_infinity'][-n_ion:])
        G = float(full[:split] @ mu_neutral + full[-n_ion:] @ mu_ions)
        return G, np.concatenate([mu_neutral, mu_ions])[self.aq_positions]

    def _lnphi(self, phase, n):
        """ln(phi) of the phase's species at amounts n (the nonideal part of mu)."""
        ns = self.ns
        if phase == 'aqueous' and self.aqueous_model == 'electrolyte':
            full = np.zeros(len(ns['AQUEOUS_SPECIES']))
            full[self.aq_positions] = n
            try:
                phi = ns['phi_electrolyte'](full)
            except ValueError as exc:
                raise EOSFailure(str(exc)) from exc
            return np.log(phi)[self.aq_positions]
        names = self.phase_species[phase]
        root = 'vapour' if phase == 'gas' else 'liquid'
        phis = np.asarray(ns['cl'].fugacity_coefficient(self._eos_model(names), ns['P'], ns['T'],
                                                        np.maximum(n, 1e-14), phase=root),
                          dtype=float).reshape(-1)
        if np.any(~np.isfinite(phis)) or np.any(phis <= 0):
            raise EOSFailure(f'Invalid PC-SAFT fugacities: phase={root}, species={list(names)}')
        return np.log(phis)

    def phase_hessian(self, phase, n):
        """d mu / d n for one fluid phase: analytic ideal part + forward-difference ln(phi) part."""
        ns = self.ns
        RT = ns['R'] * ns['T']
        if phase == 'aqueous' and self.aqueous_model == 'electrolyte':
            if n[0] <= 1e-9:
                return np.zeros((len(n), len(n)))
            # Unselected default aqueous species enter the model's mole fractions at 1e-14.
            total = np.maximum(n, 1e-14).sum() + 1e-14 * (len(ns['AQUEOUS_SPECIES']) - len(n))
        elif np.sum(n) <= 1e-13:
            return np.zeros((len(n), len(n)))
        else:
            total = np.maximum(n, 1e-14).sum()
        alive = n > 1e-14
        ideal = np.diag(np.where(alive, 1.0 / np.maximum(n, 1e-14), 0.0)) - np.outer(alive, alive) / total
        base = self._lnphi(phase, n)
        jac = np.zeros((len(n), len(n)))
        h = 1e-7 * n.sum()
        for j in range(len(n)):
            step = n.copy()
            step[j] += h
            jac[:, j] = (self._lnphi(phase, step) - base) / h
        return RT * (ideal + 0.5 * (jac + jac.T))

    def hessian_blocks(self, x):
        """Scaled Hessian of the objective as {phase: (free positions, dense block)}."""
        amounts = np.maximum(self.expand(x), 0.0)
        position = {j: k for k, j in enumerate(self.free)}
        blocks = {}
        for phase, names in self.phase_species.items():
            if not names or phase == 'solid':
                continue
            idx = self.phase_indices[phase]
            keep = [i for i, j in enumerate(idx) if j in position]
            if not keep:
                continue
            H = self.phase_hessian(phase, amounts[idx])[np.ix_(keep, keep)] / self.scale
            blocks[phase] = (np.array([position[idx[i]] for i in keep]), H)
        return blocks

    def gibbs_and_mu(self, amounts):
        """Total G [J/mol feed C] and mu for every (phase, species) pair."""
        amounts = np.maximum(amounts, 0.0)
        G, mu = 0.0, np.zeros(len(self.pairs))
        for phase, names in self.phase_species.items():
            if not names:
                continue
            idx = self.phase_indices[phase]
            n = amounts[idx]
            if phase == 'solid':
                g, m = float(n @ self.mu0[phase]), self.mu0[phase]
            elif phase == 'aqueous' and self.aqueous_model == 'electrolyte':
                g, m = self._electrolyte(n)
            else:
                g, m = self._molecular(phase, n, 'vapour' if phase == 'gas' else 'liquid')
            G += g
            mu[idx] = m
        if not np.isfinite(G) or np.any(~np.isfinite(mu)):
            raise EOSFailure('Nonfinite phase Gibbs energy.')
        return G, mu

    # ----- optimizer-facing functions -----
    def expand(self, x):
        amounts = np.zeros(len(self.pairs))
        amounts[self.free] = x
        return amounts

    def fun_and_grad(self, x, penalty=True):
        """Scaled objective and its exact gradient.

        An EOS failure returns a large penalty (penalty=True, for SLSQP whose line
        search then backtracks) or NaN (for IPOPT, which cuts the step back).
        """
        self.nfev += 1
        try:
            G, mu = self.gibbs_and_mu(self.expand(x))
        except EOSFailure as exc:
            if len(self.eos_errors) < 10:
                self.eos_errors.append(str(exc))
            if not penalty:
                return float('nan'), np.full(len(x), np.nan)
            d = x - self.x_start
            return 1e6 + float(d @ d), 2.0 * d
        self._last = (x.copy(), G, mu)
        return (G - self.G0) / self.scale, mu[self.free] / self.scale

    def check_against_model(self, x, rtol=1e-9):
        """Raise if this module's G differs from the model notebook's G at x."""
        amounts = self.expand(x)
        ours, _ = self.gibbs_and_mu(amounts)
        theirs = self.model_total_G(amounts)
        if abs(ours - theirs) > rtol * max(1.0, abs(theirs)):
            raise RuntimeError(f'htg_optimizer G={ours!r} differs from model G={theirs!r}; '
                               'the model formulas changed, update htg_optimizer')

    # ----- optimality certificate -----
    def kkt(self, x, major_tol=1e-6):
        """First-order optimality certificate with reduced costs r = (mu - A^T lambda)/RT.

        At a local minimum r_i = 0 where n_i > lower and r_i >= 0 at the bound.
        Reported:  descent_RT     = max(-r_i): a positive value is a direction that lowers G;
                   complementarity_J = RT sum (n_i - lower_i) max(r_i, 0) over minor
                                    species: G still recoverable by removing them;
                   stationarity_major_RT = max |r_i| over species above major_tol.
        """
        _, mu = self.gibbs_and_mu(self.expand(x))
        g = mu[self.free]
        RT = self.ns['R'] * self.ns['T']
        rows = [self.A_eq]
        if self.ineq_row is not None and abs(self.ineq_row @ x - self.ineq_target) < 1e-9:
            rows.append(self.ineq_row[None, :])
        A = np.vstack(rows)
        slack = x - self.lower
        major = slack > major_tol
        fit = major if major.sum() >= A.shape[0] else slack > 0
        lam, *_ = np.linalg.lstsq(A[:, fit].T, g[fit], rcond=None)
        reduced = (g - A.T @ lam) / RT
        return {'element_potentials_J': lam[:len(self.A_eq)], 'reduced_cost_RT': reduced,
                'tangent_plane_J': lam @ A,
                'descent_RT': float(max(0.0, -reduced.min())),
                'complementarity_J': float(RT * np.sum((slack * np.maximum(reduced, 0))[~major])),
                'stationarity_major_RT': float(np.max(np.abs(reduced[major]))) if major.any() else 0.0}


# ---------------------------------------------------------------- backends
def _slsqp(problem, x0, options):
    cons = [{'type': 'eq', 'fun': lambda x: problem.A_eq @ x - problem.b_eq,
             'jac': lambda x: problem.A_eq}]
    if problem.ineq_row is not None:
        cons.append({'type': 'ineq', 'fun': lambda x: problem.ineq_row @ x - problem.ineq_target,
                     'jac': lambda x: problem.ineq_row})
    bounds = [(float(lo), None) for lo in problem.lower]
    return minimize(problem.fun_and_grad, x0, jac=True, method='SLSQP', bounds=bounds,
                    constraints=cons,
                    options={'maxiter': options.get('max_iterations', 5000),
                             'ftol': options.get('ftol', 1e-12), 'disp': False})


def _trust_constr(problem, x0, options):
    cons = [LinearConstraint(problem.A_eq, problem.b_eq, problem.b_eq)]
    if problem.ineq_row is not None:
        cons.append(LinearConstraint(problem.ineq_row[None, :], problem.ineq_target, np.inf))
    return minimize(problem.fun_and_grad, x0, jac=True, method='trust-constr',
                    bounds=Bounds(problem.lower, np.inf), constraints=cons,
                    options={'maxiter': options.get('max_iterations', 5000),
                             'gtol': options.get('gtol', 1e-10), 'xtol': options.get('xtol', 1e-14),
                             'verbose': 0})


_JULIA_IPOPT = r"""
using Ipopt
function (x0, xL, xU, gL, gU, Ar, Ac, Av, Hr, Hc, pyf, pygrad, pyhess, onames, ovalues, okinds)
    x0 = Float64.(x0); xL = Float64.(xL); xU = Float64.(xU)
    gL = Float64.(gL); gU = Float64.(gU)
    Ar = Int32.(Ar); Ac = Int32.(Ac); Av = Float64.(Av); Hr = Int32.(Hr); Hc = Int32.(Hc)
    eval_f(x) = pyconvert(Float64, pyf(x))
    eval_grad_f(x, g) = (g .= pyconvert(Vector{Float64}, pygrad(x)); nothing)
    function eval_g(x, g)
        fill!(g, 0.0)
        for k in eachindex(Av)
            g[Ar[k]] += Av[k] * x[Ac[k]]
        end
        return nothing
    end
    function eval_jac_g(x, rows, cols, values)
        if values === nothing
            rows .= Ar; cols .= Ac
        else
            values .= Av
        end
        return nothing
    end
    function eval_h(x, rows, cols, sigma, lambda, values)
        if values === nothing
            rows .= Hr; cols .= Hc
        else
            values .= sigma .* pyconvert(Vector{Float64}, pyhess(x))
        end
        return nothing
    end
    prob = Ipopt.CreateIpoptProblem(length(x0), xL, xU, length(gL), gL, gU, length(Av), length(Hr),
                                    eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h)
    for (name, value, kind) in zip(onames, ovalues, okinds)
        name = pyconvert(String, name); kind = pyconvert(String, kind)
        if kind == "str"
            Ipopt.AddIpoptStrOption(prob, name, pyconvert(String, value))
        elseif kind == "int"
            Ipopt.AddIpoptIntOption(prob, name, pyconvert(Int, value))
        else
            Ipopt.AddIpoptNumOption(prob, name, pyconvert(Float64, value))
        end
    end
    prob.x .= x0
    status = Ipopt.IpoptSolve(prob)
    return (copy(prob.x), Int(status), prob.obj_val)
end
"""
_IPOPT_STATUS = {0: 'Solve_Succeeded', 1: 'Solved_To_Acceptable_Level', 2: 'Infeasible_Problem_Detected',
                 3: 'Search_Direction_Becomes_Too_Small', 4: 'Diverging_Iterates',
                 5: 'User_Requested_Stop', 6: 'Feasible_Point_Found',
                 -1: 'Maximum_Iterations_Exceeded', -2: 'Restoration_Failed',
                 -3: 'Error_In_Step_Computation', -4: 'Maximum_CpuTime_Exceeded',
                 -10: 'Not_Enough_Degrees_Of_Freedom', -11: 'Invalid_Problem_Definition',
                 -12: 'Invalid_Option', -13: 'Invalid_Number_Detected'}
_ipopt_solve = None


def _ipopt(problem, x0, options):
    """IPOPT (Ipopt.jl through juliacall) with the exact gradient and phase-block Hessian.

    An EOS failure at a trial point returns NaN, which makes IPOPT cut the step back.
    """
    from scipy.optimize import OptimizeResult
    global _ipopt_solve
    if _ipopt_solve is None:
        try:
            _ipopt_solve = problem.ns['jl'].seval(_JULIA_IPOPT)
        except Exception as exc:
            raise RuntimeError('The ipopt backend needs Ipopt.jl in the juliacall environment: '
                               "import juliapkg; juliapkg.add('Ipopt', "
                               "'b6b21f68-93f8-5de0-b562-5493be1d77c9'); juliapkg.resolve()") from exc
    n = len(x0)
    A = problem.A_eq
    gL, gU = list(problem.b_eq), list(problem.b_eq)
    if problem.ineq_row is not None:
        A = np.vstack([A, problem.ineq_row])
        gL.append(problem.ineq_target)
        gU.append(2e19)
    Ar, Ac = np.nonzero(A)
    # Lower-triangle sparsity of the phase-block Hessian (structure fixed for the solve).
    blocks = problem.hessian_blocks(x0)
    Hr, Hc = [], []
    for pos, _ in blocks.values():
        for a in range(len(pos)):
            for b in range(a + 1):
                Hr.append(max(pos[a], pos[b]))
                Hc.append(min(pos[a], pos[b]))
    nit = [0]
    last = {}

    def evaluate(x):
        x = np.array(x, dtype=float)
        if 'x' not in last or not np.array_equal(x, last['x']):
            last['x'] = x
            last['f'], last['g'] = problem.fun_and_grad(x, penalty=False)
        return last['f'], last['g']

    def f(x):
        return evaluate(x)[0]

    def grad(x):
        return evaluate(x)[1]

    def hess(x):
        nit[0] += 1
        x = np.array(x, dtype=float)
        try:
            values = []
            for pos, H in problem.hessian_blocks(x).values():
                for a in range(len(pos)):
                    for b in range(a + 1):
                        values.append(H[a, b])
            return np.array(values)
        except EOSFailure:
            return np.full(len(Hr), np.nan)

    opts = {'tol': 1e-10, 'acceptable_tol': 1e-7, 'max_iter': int(options.get('max_iterations', 3000)),
            'print_level': int(options.get('print_level', 0)),
            # The objective is O(1) and starts feasible.  A large or adaptively raised
            # barrier parameter pushes trace species into EOS-invalid compositions
            # (observed at 100-200 C with mu_strategy='adaptive').
            'mu_strategy': 'monotone', 'mu_init': 1e-6,
            'bound_relax_factor': 0.0, 'honor_original_bounds': 'yes',
            # Default bound_push (1e-2) would move trace species far from an LP vertex.
            'bound_push': 1e-6, 'bound_frac': 1e-6,
            'jac_c_constant': 'yes', 'jac_d_constant': 'yes', 'sb': 'yes'}
    opts.update(options.get('ipopt_options', {}))
    kinds = ['str' if isinstance(v, str) else 'int' if isinstance(v, int) else 'num' for v in opts.values()]
    upper = np.full(n, 2e19)
    xs, status, obj = _ipopt_solve(x0, problem.lower, upper, np.array(gL), np.array(gU),
                                   Ar + 1, Ac + 1, A[Ar, Ac], np.array(Hr, dtype=int) + 1,
                                   np.array(Hc, dtype=int) + 1, f, grad, hess,
                                   list(opts), list(opts.values()), kinds)
    status = int(status)
    return OptimizeResult(x=np.array(xs, dtype=float), fun=float(obj), status=status,
                          success=status in (0, 1), message=_IPOPT_STATUS.get(status, str(status)),
                          nit=nit[0])


BACKENDS = {'slsqp': _slsqp, 'trust-constr': _trust_constr, 'ipopt': _ipopt}


# ---------------------------------------------------------------- multistart
def _lp_point(problem, cost):
    from scipy.optimize import linprog
    extra = {}
    if problem.ineq_row is not None:
        extra = {'A_ub': -problem.ineq_row[None, :], 'b_ub': np.array([-problem.ineq_target])}
    res = linprog(cost, A_eq=problem.A_eq, b_eq=problem.b_eq,
                  bounds=[(float(lo), None) for lo in problem.lower], method='highs', **extra)
    return res.x if res.success else None


def starting_points(problem, n_random=8, seed=0):
    """Feasible starts: the given start, phase-rich/phase-lean LP vertices, random vertex mixes.

    Every point is a convex combination of LP solutions, so it satisfies A x = b,
    c x >= d and x >= lower exactly (to LP precision).
    """
    starts = [('initial', problem.x_start.copy())]
    position = {j: k for k, j in enumerate(problem.free)}
    for phase, names in problem.phase_species.items():
        if not names or phase == 'solid':
            continue
        cols = [position[j] for j in problem.phase_indices[phase] if j in position]
        if not cols:
            continue
        for sign, label in ((-1.0, 'max'), (1.0, 'min')):
            cost = np.zeros(len(problem.free))
            cost[cols] = sign
            x = _lp_point(problem, cost)
            if x is not None:
                # Blend toward the given start so no phase is exactly empty.
                starts.append((f'{label} {phase}', 0.9 * x + 0.1 * problem.x_start))
    rng = np.random.default_rng(seed)
    vertices = [x for x in (_lp_point(problem, rng.normal(size=len(problem.free)))
                            for _ in range(max(3, n_random))) if x is not None]
    for k in range(n_random if vertices else 0):
        weights = rng.dirichlet(np.ones(len(vertices)))
        starts.append((f'random {k}', np.sum([w * v for w, v in zip(weights, vertices)], axis=0)))
    return starts


# ---------------------------------------------------------------- stability
def tangent_plane_test(problem, x, tangent_plane):
    """Minimum tangent-plane distance (units of RT) for each nearly empty fluid phase.

    tpd(y) = sum_i y_i (mu_i(y) - (A^T lambda)_i) / RT over compositions y of that
    phase, with A^T lambda (``tangent_plane``, per free variable) from the KKT fit.
    A negative minimum means adding the phase would lower G, so the solution is
    not the global minimum of the selected-phase problem.  Ions are excluded from
    the aqueous trial phase, which keeps trial compositions electroneutral.
    """
    amounts = problem.expand(x)
    RT = problem.ns['R'] * problem.ns['T']
    position = {j: k for k, j in enumerate(problem.free)}
    target = np.full(len(problem.pairs), np.nan)
    target[problem.free] = tangent_plane
    report = {}
    for phase, names in problem.phase_species.items():
        if not names or phase == 'solid':
            continue
        idx = problem.phase_indices[phase]
        if amounts[idx].sum() > 1e-6:
            continue
        ions = problem.ns.get('ION_DATA', {})
        trial = [i for i, j in enumerate(idx) if j in position and names[i] not in ions]
        if not trial:
            continue

        def tpd(y_trial):
            n = np.zeros(len(idx))
            n[trial] = np.maximum(y_trial, 1e-14)
            n /= n.sum()
            try:
                if phase == 'aqueous' and problem.aqueous_model == 'electrolyte':
                    if 'water' not in names or n[names.index('water')] <= 1e-9:
                        return np.inf, np.zeros(len(trial))  # the model needs solvent water
                    _, mu = problem._electrolyte(n)
                else:
                    mu = problem._molecular_mu_at(n, names, problem.mu0[phase],
                                                  'vapour' if phase == 'gas' else 'liquid')
            except (EOSFailure, ValueError):
                return np.inf, np.zeros(len(trial))
            d = (mu - target[idx])[trial] / RT
            return float(n[trial] @ d), d
        best = np.inf
        cons = [{'type': 'eq', 'fun': lambda y: y.sum() - 1.0, 'jac': lambda y: np.ones(len(y))}]
        guesses = [np.full(len(trial), 1.0 / len(trial))] + [
            np.where(np.arange(len(trial)) == k, 1.0 - 1e-6 * (len(trial) - 1), 1e-6)
            for k in range(len(trial))]
        for y0 in guesses:
            if not np.isfinite(tpd(y0)[0]):
                continue
            res = minimize(tpd, y0, jac=True, method='SLSQP', bounds=[(0.0, 1.0)] * len(trial),
                           constraints=cons, options={'maxiter': 200, 'ftol': 1e-12})
            value = tpd(res.x)[0]
            best = min(best, value)
        report[phase] = best
    return report


def optimize(problem, options):
    """Run the selected backend (optionally multistart); return (scipy result, diagnostics).

    options: method ('slsqp' | 'trust-constr' | 'ipopt'), multistart (number of
    random starts in addition to the given and phase-rich ones; 0 = single start),
    seed, stability_test (bool), max_iterations, plus backend-specific keys.
    """
    method = options.get('method', 'ipopt')
    if method not in BACKENDS:
        raise ValueError(f"Unknown optimizer {method!r}; choose from {sorted(BACKENDS)}")
    problem.check_against_model(problem.x_start)
    t0 = time.perf_counter()
    n_random = int(options.get('multistart', 0))
    starts = (starting_points(problem, n_random, options.get('seed', 0)) if n_random > 0
              else [('initial', problem.x_start.copy())])
    runs, best, result = [], None, None
    for label, x0 in starts:
        try:
            result = BACKENDS[method](problem, x0, options)
        except Exception as exc:  # a start may lie outside the EOS domain
            if len(starts) == 1:
                raise
            runs.append({'start': label, 'G': np.nan, 'message': f'{type(exc).__name__}: {exc}'[:200]})
            continue
        G = np.nan
        if result.success and np.all(np.isfinite(result.x)):
            try:
                G, _ = problem.gibbs_and_mu(problem.expand(np.maximum(result.x, problem.lower)))
            except EOSFailure:
                pass
        runs.append({'start': label, 'G': G, 'message': str(result.message), 'nit': result.nit})
        if np.isfinite(G) and (best is None or G < best[0]):
            best = (G, result)
    if best is not None:
        result = best[1]
    elif result is None:
        raise RuntimeError(f'All optimizer starts failed: {runs}')
    # otherwise return the last (unsuccessful) result so the caller reports its message
    diagnostics = {'method': method, 'nfev': problem.nfev, 'seconds': time.perf_counter() - t0,
                   'starts': runs}
    if np.all(np.isfinite(result.x)):
        x = np.maximum(result.x, problem.lower)
        problem.check_against_model(x)
        diagnostics.update(problem.kkt(x))
        if options.get('stability_test', n_random > 0):
            diagnostics['tangent_plane_RT'] = tangent_plane_test(
                problem, x, diagnostics['tangent_plane_J'])
    return result, diagnostics


# ---------------------------------------------------------------- parallel sweeps
# juliacall cannot be used from several Python threads, so parallel solves use
# processes.  Each worker starts its own Julia and compiles Clapeyron (~40 s), so
# this pays off for sweeps of many cases, not for a single solve.
_worker_model = None


def _load_model_notebook(path):
    import json
    import types
    from pathlib import Path
    document = json.loads(Path(path).read_text(encoding='utf-8'))
    module = types.ModuleType('htg_model')
    for index, cell in enumerate(document['cells']):
        if cell['cell_type'] == 'code':
            exec(compile(''.join(cell['source']), f'{path}:cell{index + 1}', 'exec'), module.__dict__)
    return module


def _worker_init(model_notebook, database_paths, initialize_kwargs, module_overrides):
    global _worker_model
    _worker_model = _load_model_notebook(model_notebook)
    _worker_model.initialize_model(database_paths, **initialize_kwargs)
    for name, value in (module_overrides or {}).items():
        setattr(_worker_model, name, value)


def _worker_solve(settings):
    try:
        return _worker_model.solve_equilibrium(settings)
    except Exception as exc:
        return {'error': f'{type(exc).__name__}: {exc}'}


def parallel_solve(settings_list, database_paths, model_notebook='HTG_Model.ipynb', workers=None,
                   initialize_kwargs=None, module_overrides=None):
    """Solve each settings dict in a separate worker process; results keep the input order.

    module_overrides sets model attributes in every worker after initialization,
    e.g. {'AW': {...}} where a case notebook changes model.AW before solving.
    A failed case returns {'error': message} instead of raising.
    """
    import multiprocessing as mp
    import os
    from concurrent.futures import ProcessPoolExecutor
    workers = workers or min(len(settings_list), max(1, (os.cpu_count() or 2) // 2))
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn'),
                             initializer=_worker_init,
                             initargs=(str(model_notebook), {k: str(v) for k, v in database_paths.items()},
                                       initialize_kwargs or {}, module_overrides)) as pool:
        return list(pool.map(_worker_solve, settings_list))
