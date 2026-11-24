"""Results of equilibrium calculations, as objects and tables.

Amounts are in mol per mol of dry-feed carbon. Tables also give yields per kg
of original dry feed (salts are extra mass and not included in the basis).
"""
import dataclasses
import tomllib
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .constants import PHASES

PHASE_LABELS = {'gas': 'gas', 'aqueous': 'aqueous', 'organic_liquid': 'organic liquid',
                'solid': 'solid'}
CHOICES_FILE = 'choices.toml'
# Columns of alternatives.csv that are not species-table columns.
ALTERNATIVE_COLUMNS = ('result', 'gibbs_J_per_mol_C', 'n_starts', 'phase_gibbs_J_per_mol_C', 'start_labels')


@dataclass
class Alternative:
    """One distinct result at a point (the converged starts that found it).

    label         'A', 'B', ... by increasing Gibbs energy; A is the main result
    gibbs         total Gibbs energy [J/mol dry-feed C]
    amounts       {(phase, species): mol / mol dry-feed C}
    n_starts      number of starts that found it; start_labels: their names
    phase_gibbs   {phase: Gibbs energy [J/mol dry-feed C]}
    """
    label: str
    gibbs: float
    amounts: dict
    n_starts: int
    start_labels: list
    phase_gibbs: dict = field(default_factory=dict)


@dataclass
class Result:
    """Equilibrium at one temperature and pressure.

    amounts                {(phase, species): mol / mol dry-feed C}
    gibbs_J_per_mol_C      total Gibbs energy at the solution
    gas_carbon_fraction    carbon in gas / all feed carbon (including salt carbon)
    carbon_bound_violation largest violation of the carbon bounds [fraction of feed carbon]
    phase_rules            what the aqueous cut-off did at this temperature
                           ('dropped_constraints': bounds dropped at the cut-off)
    optimizer              solver diagnostics (KKT residuals, starts, results, ...)
    alternatives           the distinct results [Alternative], A first (A = this result)
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
    carbon_bound_violation: float
    phase_rules: dict
    optimizer: dict
    atomic_weights: dict
    species_data: dict = field(repr=False)
    alternatives: list = field(default_factory=list, repr=False)

    @property
    def n_results(self):
        return max(1, len(self.alternatives))

    @property
    def notes(self):
        """Short notes for the summary table, e.g. 'aqueous carbon bound dropped'."""
        return '; '.join(f'{name} bound dropped'
                         for name in (self.phase_rules or {}).get('dropped_constraints', []))

    def for_alternative(self, label):
        """This result with the amounts and Gibbs energies of alternative ``label``."""
        if label == 'A' and not any(a.label == 'A' for a in self.alternatives):
            return self
        alt = next((a for a in self.alternatives if a.label == label), None)
        if alt is None:
            raise KeyError(f'No result {label!r} at {self.temperature_C:g} °C, {self.pressure_bar:g} bar.')
        amounts = {pair: float(alt.amounts.get(pair, 0.0)) for pair in self.amounts}
        feed_C = self.feed_atoms[self.elements.index('C')]
        gas_C = sum(v * self.species_data[pair]['elements'].get('C', 0) for pair, v in amounts.items()
                    if pair[0] == 'gas')
        product = np.array([sum(v * self.species_data[pair]['elements'].get(e, 0)
                                for pair, v in amounts.items()) for e in self.elements], dtype=float)
        phase_gibbs = {phase: 0.0 for phase in self.phase_gibbs_J_per_mol_C}
        phase_gibbs.update({k: float(v) for k, v in alt.phase_gibbs.items()})
        return dataclasses.replace(self, amounts=amounts, gibbs_J_per_mol_C=float(alt.gibbs),
                                   phase_gibbs_J_per_mol_C=phase_gibbs, product_atoms=product,
                                   gas_carbon_fraction=float(gas_C / feed_C))

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
                'descent_RT': self.optimizer.get('descent_RT', np.nan),
                'n_results': self.n_results, 'carbon_bound_violation': self.carbon_bound_violation,
                'notes': self.notes, 'error': ''}

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
    user_species: dict | None = field(default=None, repr=False)  # see htg.user_species.run_record

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

    def alternatives_table(self):
        """Long table of the distinct results at points with more than one result.

        Columns: T_C, P_bar, result, gibbs_J_per_mol_C, n_starts, phase, species,
        mol_per_mol_C, then the other species-table columns, the Gibbs energy of
        the row's phase and the starts that found the result.
        """
        first = ['T_C', 'P_bar', 'result', 'gibbs_J_per_mol_C', 'n_starts', 'phase', 'species',
                 'mol_per_mol_C']
        tables = []
        for r in self.results:
            if len(r.alternatives) < 2:
                continue
            for alt in r.alternatives:
                table = r.for_alternative(alt.label).species_table()
                table.insert(2, 'result', alt.label)
                table.insert(3, 'gibbs_J_per_mol_C', float(alt.gibbs))
                table.insert(4, 'n_starts', int(alt.n_starts))
                table['phase_gibbs_J_per_mol_C'] = [float(alt.phase_gibbs.get(p, np.nan)) for p in table.phase]
                table['start_labels'] = '; '.join(alt.start_labels)
                tables.append(table)
        if not tables:
            return pd.DataFrame(columns=first)
        table = pd.concat(tables, ignore_index=True)
        return table[first + [c for c in table.columns if c not in first]]

    def save(self, directory):
        """Write summary.csv, species.csv, phases.csv, alternatives.csv and case.toml
        (plus the comparison data, if the case has any) to ``directory``."""
        import shutil
        from pathlib import Path
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        self.summary().to_csv(directory / 'summary.csv', index=False)
        self.species_table().to_csv(directory / 'species.csv', index=False)
        self.phase_summary().to_csv(directory / 'phases.csv', index=False)
        self.alternatives_table().to_csv(directory / 'alternatives.csv', index=False)
        saved = self.case.replace()
        if self.case.comparison and Path(self.case.comparison['file']).exists():
            shutil.copyfile(self.case.comparison['file'], directory / 'comparison.csv')
            saved.comparison = {**self.case.comparison, 'file': 'comparison.csv'}
        saved.save(directory / 'case.toml')
        if self.user_species is not None:
            from .user_species import write_run_record
            write_run_record(directory, self.user_species)
        return directory


# ---------------------------------------------------------------- chosen results
def _point(T, P):
    return float(T), float(P)


def read_choices(folder):
    """{(T_C, P_bar): result label} from choices.toml in ``folder`` ({} without the file).

    The file holds only points that do not use result A:
        [[point]]
        T_C = 400
        P_bar = 240
        result = "B"
    """
    path = Path(folder) / CHOICES_FILE
    if not path.exists():
        return {}
    data = tomllib.loads(path.read_text())
    return {_point(entry['T_C'], entry['P_bar']): str(entry['result'])
            for entry in data.get('point', []) if str(entry.get('result', 'A')) != 'A'}


def write_choices(folder, choices):
    """Write choices.toml (points that use A are left out; no file when every point uses A)."""
    path = Path(folder) / CHOICES_FILE
    chosen = {_point(*point): str(label) for point, label in choices.items() if str(label) != 'A'}
    if not chosen:
        path.unlink(missing_ok=True)
        return None
    lines = ['# The result each point uses in plots and tables (points not listed use A).', '']
    for (T, P), label in sorted(chosen.items()):
        lines += ['[[point]]', f'T_C = {T!r}', f'P_bar = {P!r}', f'result = "{label}"', '']
    path.write_text('\n'.join(lines))
    return path


def _at(table, T, P):
    mask = np.isclose(table.T_C.astype(float), T)
    if 'P_bar' in table:
        mask &= np.isclose(table.P_bar.astype(float), P)
    return mask


def apply_choices(summary, species, phases, alternatives, choices):
    """Tables with the chosen result at each point: (summary, species, phases).

    alternatives  the alternatives table (alternatives.csv); choices {(T, P): label}.
    The summary gains a 'result' column when there are alternatives. At a point
    that uses another result than A, its amounts replace the species rows, the
    phase rows are recomputed from them, and the summary's Gibbs energy and gas
    carbon follow; the solver columns stay those of A. Unknown labels are ignored.
    """
    summary, species, phases = summary.copy(), species.copy(), phases.copy()
    if alternatives is None or alternatives.empty:
        return summary, species, phases
    if 'result' not in summary:
        status = summary['status'] if 'status' in summary else pd.Series('solved', index=summary.index)
        summary.insert(min(3, len(summary.columns)), 'result', np.where(status == 'solved', 'A', ''))
    species_columns = [c for c in alternatives.columns if c not in ALTERNATIVE_COLUMNS and c in species]
    for point, label in choices.items():
        T, P = _point(*point)
        rows = alternatives[_at(alternatives, T, P) & (alternatives.result == label)]
        if rows.empty or label == 'A':
            continue
        mask = _at(species, T, P)
        if mask.sum() != len(rows):
            continue
        species.loc[mask, species_columns] = rows[species_columns].to_numpy()
        carbon = rows.groupby('phase').carbon_mol_per_kg_dry_feed.sum()
        total_carbon = float(carbon.sum())
        gibbs = rows.groupby('phase').phase_gibbs_J_per_mol_C.first()
        for index in phases.index[_at(phases, T, P)]:
            phase = phases.at[index, 'phase']
            part = rows[rows.phase == phase]
            phases.at[index, 'mol_per_mol_C'] = float(part.mol_per_mol_C.sum())
            if 'g_per_kg_dry_feed' in phases:
                phases.at[index, 'g_per_kg_dry_feed'] = float(part.g_per_kg_dry_feed.sum())
            phases.at[index, 'carbon_percent'] = (100 * float(carbon.get(phase, 0.0)) / total_carbon
                                                  if total_carbon else 0.0)
            if 'gibbs_J_per_mol_C' in phases:
                value = gibbs.get(phase, np.nan)
                phases.at[index, 'gibbs_J_per_mol_C'] = 0.0 if pd.isna(value) else float(value)
        row = _at(summary, T, P)
        summary.loc[row, 'result'] = label
        summary.loc[row, 'gibbs_J_per_mol_C'] = float(rows.gibbs_J_per_mol_C.iloc[0])
        if 'gas_carbon_percent' in summary:
            summary.loc[row, 'gas_carbon_percent'] = (100 * float(carbon.get('gas', 0.0)) / total_carbon
                                                      if total_carbon else 0.0)
    return summary, species, phases
