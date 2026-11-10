"""The Model: databases + species registry, and the calculations that use them."""
import os
from pathlib import Path

import pandas as pd

from . import equilibrium
from .databases import Databases
from .species import Registry

REPO_DATA = Path(__file__).resolve().parents[2] / 'data'


def default_data_dir():
    """$HTG_DATA_DIR if set, else the repository's data/ folder."""
    return Path(os.environ.get('HTG_DATA_DIR', REPO_DATA))


class Model:
    """Thermodynamic model built from the databases in ``data_dir``.

    >>> model = htg.Model.load()
    >>> result = model.solve(case)          # one temperature
    >>> sweep = model.sweep(case)           # all temperatures of the case
    >>> model.species()                     # species available for case files
    """

    def __init__(self, databases, registry, data_dir=None):
        self.databases = databases
        self.registry = registry
        self.data_dir = data_dir

    @classmethod
    def load(cls, data_dir=None, species_file=None, nasa_names=None, reference_temperature_C=350.0):
        """Read the databases and discover species.

        species_file  registry file (default: data_dir/species.toml)
        nasa_names    replaces the [nasa_names] table of the species file
        reference_temperature_C  used only to evaluate temperature-dependent SAFT
                      parameters while checking which species have complete data
        """
        data_dir = Path(data_dir) if data_dir is not None else default_data_dir()
        databases = Databases.load(data_dir)
        registry = Registry.from_file(species_file or data_dir / 'species.toml', nasa_names)
        registry.discover(databases, float(reference_temperature_C) + 273.15)
        return cls(databases, registry, data_dir)

    # ---------------------------------------------------------------- calculations
    def solve(self, case, temperature_C=None):
        """Equilibrium for ``case`` at one temperature (the case's, or ``temperature_C``)."""
        if temperature_C is not None:
            case = case.at(temperature_C)
        return equilibrium.solve(self, case)

    def sweep(self, case, progress=print, workers=1):
        """Equilibrium at every temperature of the case; see htg.sweep.run_sweep."""
        from .sweep import run_sweep
        return run_sweep(self, case, progress=progress, workers=workers)

    def check(self, case):
        """Preflight every temperature without running the equation of state.

        Checks species names and data coverage, temperature ranges of the
        standard-state data and feasibility of all balances and constraints.
        """
        from .sweep import preflight
        return preflight(self, case)

    # ---------------------------------------------------------------- species information
    def species(self, include_unavailable=False, phase=None, search=None):
        """Species catalogue (ready species only, unless include_unavailable)."""
        table = self.registry.available(include_unavailable)
        if phase:
            table = table[table.Phase.str.contains(phase, case=False)]
        if search:
            text = str(search).casefold()
            table = table[table.Species.astype(str).str.casefold().str.contains(text, regex=False)
                          | table.Formula.astype(str).str.casefold().str.contains(text, regex=False)]
        return table.reset_index(drop=True)

    def mineral_inventory(self, temperature_C=None, elements=None):
        """SUPCRT minerals available as pure solids, with hydrate flag and T limit."""
        rows = []
        for name, d in self.registry.solid.items():
            if d.get('adapter') != 'supcrt_mineral':
                continue
            if elements and not set(elements).intersection(d['elements']):
                continue
            high = d['supcrt_mineral']['Tmax']
            rows.append({'Species': name, 'Formula': d['formula'], 'Elements': dict(d['elements']),
                         'Hydrate': d['hydrate'], 'Tmax_K': high,
                         'Within_T_range': temperature_C is None or 298.15 <= temperature_C + 273.15 <= high,
                         'References_in_XML': d['references'], 'Date_in_XML': d['date']})
        return pd.DataFrame(rows)

    def product_candidates(self, elements, temperature_C=None):
        """Neutral gas molecules and pure solids that can carry ``elements`` (e.g. salt elements)."""
        registry, nasa_db = self.registry, self.databases.nasa
        target = set(elements)
        rows = []
        for phase, entries in [('gas', registry.all_molecular), ('solid', registry.solid)]:
            for name, d in entries.items():
                if name in registry.ion or not target.intersection(d['elements']):
                    continue
                reason = ''
                if d.get('adapter') == 'supcrt_mineral':
                    high = d['supcrt_mineral']['Tmax']
                    if temperature_C is not None and not 298.15 <= temperature_C + 273.15 <= high:
                        reason = 'Outside SUPCRT mineral temperature range'
                    source = 'supcrt92.xml'
                else:
                    record = nasa_db.get(d.get('nasa'), {})
                    if record.get('phase') != ('S' if phase == 'solid' else 'G'):
                        continue
                    source = 'NASA-7 database'
                    if temperature_C is not None and not (
                            record['T_low'] <= temperature_C + 273.15 <= record['T_high']):
                        reason = 'Outside NASA temperature range'
                if phase == 'gas':
                    T = (temperature_C if temperature_C is not None else 350.0) + 273.15
                    try:
                        registry.resolve_parameters(self.databases, T, [name])
                    except (ValueError, KeyError) as exc:
                        reason = str(exc)
                rows.append({'Species': name, 'Phase': phase, 'Elements': dict(d['elements']),
                             'Source': source, 'Ready': not reason, 'Reason': reason})
        return pd.DataFrame(rows, columns=['Species', 'Phase', 'Elements', 'Source', 'Ready', 'Reason'])
