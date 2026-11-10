"""Results of equilibrium calculations, as objects and tables.

Amounts are in mol per mol of dry-feed carbon. Tables also give yields per kg
of original dry feed (salts are extra mass and not included in the basis).
"""
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .constants import PHASES

PHASE_LABELS = {'gas': 'gas', 'aqueous': 'aqueous', 'organic_liquid': 'organic liquid',
                'solid': 'solid'}


@dataclass
class Result:
    """Equilibrium at one temperature and pressure.

    amounts                {(phase, species): mol / mol dry-feed C}
    gibbs_J_per_mol_C      total Gibbs energy at the solution
    gas_carbon_fraction    carbon in gas / all feed carbon (including salt carbon)
    phase_rules            what the aqueous cut-off did at this temperature
    optimizer              solver diagnostics (KKT residuals, multistart runs, ...)
    """
    temperature_C: float
    pressure_bar: float
    amounts: dict
    initial_amounts: dict
    gibbs_J_per_mol_C: float
    phase_gibbs_J_per_mol_C: dict
    balance_residual: float
    charge_residual: float
    iterations: int
    eos_trial_errors: list
    feed_C_mol_per_g_dry: float
    elements: list
    feed_atoms: np.ndarray
    product_atoms: np.ndarray
    gas_carbon_fraction: float
    gas_carbon_min_violation: float
    phase_rules: dict
    optimizer: dict
    atomic_weights: dict
    species_data: dict = field(repr=False)

    def species_table(self):
        """One row per (phase, species): amounts, yields per kg dry feed, carbon, molality."""
        factor = self.feed_C_mol_per_g_dry * 1000.0  # mol/mol C -> mol/kg dry feed
        aw = self.atomic_weights
        water_kg = self.amounts.get(('aqueous', 'water'), 0.0) * (2 * aw['H'] + aw['O']) / 1000.0
        rows = []
        for (phase, species), amount in self.amounts.items():
            data = self.species_data[(phase, species)]
            is_ion = bool(data['charge'])
            mw = sum(aw[e] * c for e, c in data['elements'].items())
            rows.append({
                'T_C': self.temperature_C, 'P_bar': self.pressure_bar, 'phase': phase,
                'group': 'aqueous ions' if is_ion else 'aqueous neutrals' if phase == 'aqueous' else phase,
                'species': species, 'formula': data['formula'],
                'mol_per_mol_C': amount, 'mol_per_kg_dry_feed': amount * factor,
                'g_per_kg_dry_feed': amount * factor * mw,
                'carbon_mol_per_kg_dry_feed': amount * factor * data['elements'].get('C', 0.0),
                'molality_mol_per_kg_water': amount / water_kg if is_ion and water_kg > 0 else np.nan,
            })
        return pd.DataFrame(rows)

    def phase_summary(self):
        """Per phase: total amount, mass, carbon share of all feed carbon, Gibbs energy."""
        table = self.species_table()
        feed_C = self.feed_atoms[self.elements.index('C')]
        rows = []
        for phase in PHASES:
            part = table[table.phase == phase]
            carbon = float((part.mol_per_mol_C * [self.species_data[(phase, s)]['elements'].get('C', 0)
                                                  for s in part.species]).sum()) if len(part) else 0.0
            rows.append({'T_C': self.temperature_C, 'P_bar': self.pressure_bar, 'phase': phase,
                         'mol_per_mol_C': float(part.mol_per_mol_C.sum()),
                         'g_per_kg_dry_feed': float(part.g_per_kg_dry_feed.sum()),
                         'carbon_percent': 100 * carbon / feed_C,
                         'gibbs_J_per_mol_C': self.phase_gibbs_J_per_mol_C.get(phase, 0.0)})
        return pd.DataFrame(rows)

    def summary(self):
        return {'T_C': self.temperature_C, 'P_bar': self.pressure_bar, 'status': 'solved',
                'gibbs_J_per_mol_C': self.gibbs_J_per_mol_C,
                'gas_carbon_percent': 100 * self.gas_carbon_fraction,
                'balance_residual': self.balance_residual, 'charge_residual': self.charge_residual,
                'iterations': self.iterations,
                'aqueous_removed': self.phase_rules.get('aqueous_disabled'),
                'solver': self.optimizer.get('method'),
                'descent_RT': self.optimizer.get('descent_RT', np.nan), 'error': ''}

    def __repr__(self):
        major = sorted(((v, p, s) for (p, s), v in self.amounts.items() if v > 1e-3), reverse=True)
        lines = [f'Result at {self.temperature_C:g} °C, {self.pressure_bar:g} bar: '
                 f'G = {self.gibbs_J_per_mol_C:.6g} J/mol C, '
                 f'{100 * self.gas_carbon_fraction:.1f} % of carbon in gas']
        lines += [f'  {p:15s} {s:20s} {v:.4g} mol/mol C' for v, p, s in major[:12]]
        return '\n'.join(lines)


@dataclass
class SweepResult:
    """Results over a temperature-pressure grid.

    points  [(T_C, P_bar, Result | error message)] in grid order
    """
    case: object
    points: list

    @property
    def results(self):
        return [r for *_, r in self.points if isinstance(r, Result)]

    def summary(self):
        rows = []
        for T, P, r in self.points:
            if isinstance(r, Result):
                rows.append(r.summary())
            else:
                rows.append({'T_C': T, 'P_bar': P, 'status': 'failed', 'error': r})
        return pd.DataFrame(rows)

    def species_table(self):
        tables = [r.species_table() for r in self.results]
        return pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()

    def phase_summary(self):
        tables = [r.phase_summary() for r in self.results]
        return pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()

    def save(self, directory):
        """Write summary.csv, species.csv, phases.csv and case.toml (plus the
        comparison data, if the case has any) to ``directory``."""
        import shutil
        from pathlib import Path
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        self.summary().to_csv(directory / 'summary.csv', index=False)
        self.species_table().to_csv(directory / 'species.csv', index=False)
        self.phase_summary().to_csv(directory / 'phases.csv', index=False)
        saved = self.case.replace()
        if self.case.comparison and Path(self.case.comparison['file']).exists():
            shutil.copyfile(self.case.comparison['file'], directory / 'comparison.csv')
            saved.comparison = {**self.case.comparison, 'file': 'comparison.csv'}
        saved.save(directory / 'case.toml')
        return directory
