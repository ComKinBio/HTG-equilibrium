"""Optimizer backends, multistart and optimality diagnostics for the Gibbs minimization.

The phase models give G and its exact gradient: for every phase G = sum n_i mu_i
with mu_i from the same Helmholtz model, so dG/dn_i = mu_i (no finite
differences). IPOPT also receives a phase-block Hessian.

    method        'ipopt' (default; Ipopt.jl via juliacall), 'slsqp', 'trust-constr'
    multistart    k > 0 adds phase-rich LP starts and k random feasible starts;
                  the lowest G wins
    stability_test tangent-plane test of absent fluid phases (default: on with multistart)
"""
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, OptimizeResult, linprog, minimize

from .backends import julia
from .constants import R
from .phases import EOSFailure


class GibbsProblem:
    """min_x (G(expand(x)) - G0)/scale  s.t.  A x = b,  c x >= d,  x >= lower.

    ``x`` holds the free (phase, species) amounts.
    """

    def __init__(self, phases, pairs, free, phase_indices, ions, T, A_eq, b_eq, ineq_row,
                 ineq_target, lower, x_start, G0, scale, model_total_G):
        self.phases = phases
        self.pairs = list(pairs)
        self.free = np.asarray(free, dtype=int)
        self.phase_indices = phase_indices
        self.ions = ions
        self.T = float(T)
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

    def fluid_phases(self):
        return [phase for phase in self.phases if phase != 'solid']

    def expand(self, x):
        amounts = np.zeros(len(self.pairs))
        amounts[self.free] = x
        return amounts

    def gibbs_and_mu(self, amounts):
        """Total G [J/mol feed C] and mu for every (phase, species) pair."""
        amounts = np.maximum(amounts, 0.0)
        G, mu = 0.0, np.zeros(len(self.pairs))
        for phase, model in self.phases.items():
            idx = self.phase_indices[phase]
            g, m = model.gibbs_and_mu(amounts[idx])
            G += g
            mu[idx] = m
        if not np.isfinite(G) or np.any(~np.isfinite(mu)):
            raise EOSFailure('Nonfinite phase Gibbs energy.')
        return G, mu

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
        return (G - self.G0) / self.scale, mu[self.free] / self.scale

    def hessian_blocks(self, x):
        """Scaled Hessian of the objective as {phase: (free positions, dense block)}."""
        amounts = np.maximum(self.expand(x), 0.0)
        position = {j: k for k, j in enumerate(self.free)}
        blocks = {}
        for phase in self.fluid_phases():
            idx = self.phase_indices[phase]
            keep = [i for i, j in enumerate(idx) if j in position]
            if not keep:
                continue
            H = self.phases[phase].hessian(amounts[idx])[np.ix_(keep, keep)] / self.scale
            blocks[phase] = (np.array([position[idx[i]] for i in keep]), H)
        return blocks

    def check_against_model(self, x, rtol=1e-9):
        """Raise if the optimizer's G (gibbs_and_mu) differs from the model's G (gibbs)."""
        amounts = self.expand(x)
        ours, _ = self.gibbs_and_mu(amounts)
        theirs = self.model_total_G(amounts)
        if abs(ours - theirs) > rtol * max(1.0, abs(theirs)):
            raise RuntimeError(f'Optimizer G={ours!r} differs from model G={theirs!r}')

    def kkt(self, x, major_tol=1e-6):
        """First-order optimality certificate with reduced costs r = (mu - A^T lambda)/RT.

        At a local minimum r_i = 0 where n_i > lower and r_i >= 0 at the bound.
            descent_RT            max(-r_i): positive means a direction that lowers G
            complementarity_J     RT sum (n_i - lower_i) max(r_i, 0) over minor species
            stationarity_major_RT max |r_i| over species above major_tol
        """
        _, mu = self.gibbs_and_mu(self.expand(x))
        g = mu[self.free]
        RT = R * self.T
        rows = [self.A_eq]
        # Interior-point solutions end slightly inside an active bound, hence 1e-6.
        if self.ineq_row is not None and abs(self.ineq_row @ x - self.ineq_target) < 1e-6:
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
    """IPOPT with the exact gradient and phase-block Hessian.

    An EOS failure at a trial point returns NaN, which makes IPOPT cut the step back.
    """
    global _ipopt_solve
    if _ipopt_solve is None:
        try:
            _ipopt_solve = julia().seval(_JULIA_IPOPT)
        except Exception as exc:
            raise RuntimeError('The ipopt solver needs Ipopt.jl in the Julia environment. It is '
                               'declared in src/htg/juliapkg.json; reinstall with `uv sync` and '
                               'restart Python.') from exc
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
            # The objective is O(1) and starts feasible. A large or adaptively raised
            # barrier parameter pushes trace species into EOS-invalid compositions
            # (observed at 100-200 C with mu_strategy='adaptive'). mu_init = 1e-6
            # sometimes stopped at non-stationary points (manure with ions and
            # minerals, 100-360 C); 1e-9 reached the lowest G found by any setting
            # at every point of the reference cases and example sweeps.
            'mu_strategy': 'monotone', 'mu_init': 1e-9,
            'bound_relax_factor': 0.0, 'honor_original_bounds': 'yes',
            # Default bound_push (1e-2) would move trace species far from an LP vertex.
            'bound_push': 1e-6, 'bound_frac': 1e-6,
            'jac_c_constant': 'yes', 'jac_d_constant': 'yes', 'sb': 'yes'}
    opts.update(options.get('ipopt_options', {}))
    kinds = ['str' if isinstance(v, str) else 'int' if isinstance(v, int) else 'num' for v in opts.values()]
    xs, status, obj = _ipopt_solve(x0, problem.lower, np.full(n, 2e19), np.array(gL), np.array(gU),
                                   Ar + 1, Ac + 1, A[Ar, Ac], np.array(Hr, dtype=int) + 1,
                                   np.array(Hc, dtype=int) + 1, lambda x: evaluate(x)[0],
                                   lambda x: evaluate(x)[1], hess, list(opts), list(opts.values()),
                                   kinds)
    status = int(status)
    return OptimizeResult(x=np.array(xs, dtype=float), fun=float(obj), status=status,
                          success=status in (0, 1), message=_IPOPT_STATUS.get(status, str(status)),
                          nit=nit[0])


BACKENDS = {'slsqp': _slsqp, 'trust-constr': _trust_constr, 'ipopt': _ipopt}


# ---------------------------------------------------------------- multistart
def _lp_point(problem, cost):
    extra = {}
    if problem.ineq_row is not None:
        extra = {'A_ub': -problem.ineq_row[None, :], 'b_ub': np.array([-problem.ineq_target])}
    res = linprog(cost, A_eq=problem.A_eq, b_eq=problem.b_eq,
                  bounds=[(float(lo), None) for lo in problem.lower], method='highs', **extra)
    return res.x if res.success else None


def starting_points(problem, n_random=8, seed=0):
    """Feasible starts: the given start, phase-rich/phase-lean LP vertices, random vertex mixes.

    Every point is a convex combination of LP solutions, so it satisfies the
    constraints exactly (to LP precision).
    """
    starts = [('initial', problem.x_start.copy())]
    position = {j: k for k, j in enumerate(problem.free)}
    for phase in problem.fluid_phases():
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
    phase. A negative minimum means adding the phase would lower G, so the
    solution is not the global minimum of the selected-phase problem. Ions are
    excluded from the aqueous trial phase (keeps trial compositions electroneutral).
    """
    amounts = problem.expand(x)
    RT = R * problem.T
    position = {j: k for k, j in enumerate(problem.free)}
    target = np.full(len(problem.pairs), np.nan)
    target[problem.free] = tangent_plane
    report = {}
    for phase in problem.fluid_phases():
        model = problem.phases[phase]
        names = model.species
        idx = problem.phase_indices[phase]
        if amounts[idx].sum() > 1e-6:
            continue
        trial = [i for i, j in enumerate(idx) if j in position and names[i] not in problem.ions]
        if not trial:
            continue
        electrolyte = hasattr(model, 'all_species')

        def tpd(y_trial):
            n = np.zeros(len(idx))
            n[trial] = np.maximum(y_trial, 1e-14)
            n /= n.sum()
            try:
                if electrolyte and ('water' not in names or n[names.index('water')] <= 1e-9):
                    return np.inf, np.zeros(len(trial))  # the model needs solvent water
                mu = model.mu_at(n)
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
            best = min(best, tpd(res.x)[0])
        report[phase] = best
    return report


POLISH_THRESHOLD_RT = 0.1


def _polish(problem, result, G, options):
    """Check IPOPT's answer; polish it with SLSQP when it is not stationary.

    IPOPT can report success at a point where species with clearly positive
    amounts still have reduced costs of several RT, i.e. not a local minimum
    (observed for the ions-and-minerals manure case at 250 C with
    mu_init = 1e-6). SLSQP with exact gradients, started from that point, then
    reaches a lower G.
    The lower of the two results is kept.
    """
    x = np.maximum(result.x, problem.lower)
    threshold = float(options.get('polish_threshold_RT', POLISH_THRESHOLD_RT))
    stationarity = problem.kkt(x)['stationarity_major_RT']
    info = {'stationarity_major_RT': stationarity, 'threshold_RT': threshold, 'used': False}
    if stationarity <= threshold:
        return result, info
    try:
        polished = _slsqp(problem, x, options)
        G_new = (problem.gibbs_and_mu(problem.expand(np.maximum(polished.x, problem.lower)))[0]
                 if polished.success and np.all(np.isfinite(polished.x)) else np.nan)
    except Exception as exc:
        info['message'] = f'{type(exc).__name__}: {exc}'[:200]
        return result, info
    info.update(G_ipopt=float(G), G_polished=float(G_new), message=str(polished.message))
    if np.isfinite(G_new) and G_new < G:
        info['used'] = True
        polished.nit = result.nit + polished.nit
        return polished, info
    return result, info


def optimize(problem, options):
    """Run the selected backend (optionally multistart); return (scipy result, diagnostics)."""
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
            result = None
            failure = exc
        if method == 'ipopt' and options.get('fallback', True) and (result is None or not result.success):
            # IPOPT can stop with a step-computation error when many trial points
            # fail in the equation of state; SLSQP from the same start often succeeds.
            message = str(result.message) if result is not None else f'{type(failure).__name__}: {failure}'
            try:
                retry = _slsqp(problem, x0, options)
                if retry.success:
                    result, label = retry, f'{label} (ipopt failed: {message[:80]}; slsqp)'
            except Exception:
                pass
        if result is None:
            if len(starts) == 1:
                raise failure
            runs.append({'start': label, 'G': np.nan,
                         'message': f'{type(failure).__name__}: {failure}'[:200]})
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
    if method == 'ipopt' and best is not None and options.get('polish', True):
        result, diagnostics['polish'] = _polish(problem, result, best[0], options)
    if np.all(np.isfinite(result.x)):
        x = np.maximum(result.x, problem.lower)
        problem.check_against_model(x)
        diagnostics.update(problem.kkt(x))
        if options.get('stability_test', n_random > 0):
            diagnostics['tangent_plane_RT'] = tangent_plane_test(problem, x, diagnostics['tangent_plane_J'])
    return result, diagnostics
