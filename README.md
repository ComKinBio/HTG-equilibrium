# HTG equilibrium

Constrained Gibbs-energy minimization for hydrothermal gasification (HTG):
PC-SAFT / ePC-SAFT fluid phases (via [Clapeyron.jl](https://github.com/ClapeyronThermo/Clapeyron.jl)),
NASA-7 standard states, SUPCRT92/HKF aqueous ions and minerals.

## Setup

1. Install [uv](https://docs.astral.sh/uv/getting-started/installation/).
2. From the repository root:

   ```sh
   uv sync                          # creates .venv with the locked Python packages
   uv run nbstripout --install      # strip notebook outputs on commit (once per clone)
   uv run jupyter lab               # or open the folder in VS Code and pick the .venv kernel
   ```

The first model initialization downloads Julia (if it is not already on your
`PATH`) and installs and precompiles the Julia packages into `.venv/julia_env`.
This takes several minutes once; later runs reuse it.

## Contents

| File | Purpose |
|---|---|
| `HTG_Model.ipynb` | Thermodynamic model and constrained Gibbs solver (loaded by the case notebooks) |
| `htg_optimizer.py` | Optional optimizer backends for `solve_equilibrium` (IPOPT, exact-gradient SLSQP/trust-constr, multistart, stability test), selected with `settings['optimizer']` |
| `HTG_Cases.ipynb` | Cellulose and Yakaboylu et al. (2013) Figures 7–8 cases |
| `HTG_Configurable_Case.ipynb` | Standalone configurable temperature-sweep case |
| `Thermdat_WasteMech.txt` | NASA-7 thermochemistry |
| `PROJECT_PCSAFT.csv`, `PROJECT_EPCSAFT.csv` | PC-SAFT and ePC-SAFT parameters |
| `supcrt92.xml` | SUPCRT92 aqueous species (HKF) and minerals |
| `Cellulose.xlsx` | Cellulose feed and earlier MATLAB results |
| `juliapkg.json` | Extra Julia packages (Ipopt.jl); read by juliacall when Python runs from the repository root |

Notebooks are committed without outputs; the `nbstripout` git filter removes
them automatically.
