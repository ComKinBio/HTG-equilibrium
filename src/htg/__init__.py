"""HTG equilibrium: constrained Gibbs-energy minimization for hydrothermal gasification.

    import htg
    model = htg.Model.load()                       # databases in data/
    case = htg.load_case('cases/fig78.toml')       # or htg.Case(...)
    result = model.solve(case, temperature_C=450)  # one temperature
    sweep = model.sweep(case)                      # every temperature of the case
"""
from .case import Case, Feed, Salt, load_case
from .equilibrium import SolveError
from .model import Model, default_data_dir
from .results import Result, SweepResult

__all__ = ['Case', 'Feed', 'Salt', 'load_case', 'Model', 'default_data_dir', 'Result',
           'SweepResult', 'SolveError']
__version__ = '0.2.0'
