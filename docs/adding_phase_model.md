# Adding a phase model

A *phase model* computes the Gibbs energy of one phase from its composition.
The existing models are in `src/htg/phases.py`:

| Class | Used for | Model |
|---|---|---|
| `FluidPhase` | gas, organic liquid, molecular aqueous | PC-SAFT fugacity coefficients |
| `ElectrolytePhase` | aqueous with `model = "electrolyte"` | ePC-SAFT + HKF ion standard states |
| `SolidPhase` | solids | pure phases, G = Σ nᵢ g°ᵢ |

Adding a model (for example a different equation of state for the aqueous
phase) needs some Python, but only in two places.

## 1. Write the class

A phase model is an object with these attributes and methods. Amounts `n` are
a NumPy array in the order of `species`, in mol per mol of dry-feed carbon.
Energies are in J (per mol feed carbon).

```python
class MyPhase:
    species: list[str]          # names of the phase's species, in order

    def gibbs(self, n):
        """Total Gibbs energy G = Σ nᵢ μᵢ of the phase."""

    def gibbs_and_mu(self, n):
        """(G, μ): Gibbs energy and the chemical potentials μᵢ = ∂G/∂nᵢ."""

    def mu_at(self, n):
        """Chemical potentials only (used by the stability test)."""

    def hessian(self, n):
        """∂²G/∂nᵢ∂nⱼ (used by IPOPT). Optional for the SLSQP solvers."""
```

Rules that the optimizer relies on:

- μ must be the exact derivative of G. `tests` should check this with finite
  differences for your class.
- G must be homogeneous of degree one (doubling every nᵢ doubles G).
- An empty phase (Σn ≈ 0) has G = 0; `gibbs_and_mu` should then return the
  pure-component chemical potentials (see `FluidPhase.gibbs_and_mu`).
- When the model cannot be evaluated at a composition (for example no valid
  density root), raise `htg.phases.EOSFailure` with a clear message. The
  optimizer treats this as a rejected trial point.

`FluidPhase` is the simplest complete example: μᵢ = g°ᵢ + RT ln(xᵢ φᵢ P/P°),
with the fugacity coefficients φᵢ from Clapeyron.jl. A new fluid model can
often reuse it by passing a different Clapeyron model as `eos_model`.

## 2. Build it in `equilibrium.setup`

`setup()` in `src/htg/equilibrium.py` creates one model per phase (search for
`phase_models`). Add a branch that builds your class, selected by a case
setting. For the aqueous phase, the natural switch is the existing
`aqueous_model` option (`[aqueous] model = ...` in case files):

1. allow the new value in `Case.validate` (`src/htg/case.py`);
2. build your model in `setup()` when `case.aqueous_model` has that value.

Keep the existing models unchanged and make the new one opt-in, so the
reference cases keep their results.

## 3. Test it

- Unit tests for the class: μ against finite differences of G,
  homogeneity, and known values from the literature or another program.
- A case file that uses it, checked with `uv run htg run`.
- `uv run pytest`: the reference tests must still pass, which shows that the
  existing models were not affected.
