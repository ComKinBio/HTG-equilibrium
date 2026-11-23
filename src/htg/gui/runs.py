"""Result tables of runs, live (from the worker pool) or saved (results folders)."""
from dataclasses import dataclass, field, replace
from datetime import datetime
from pathlib import Path

import pandas as pd

from ..case import load_case
from ..plotting import read_csv
from ..results import apply_choices, read_choices


@dataclass
class RunTables:
    """Everything the result and comparison views need from one run.

    summary, species and phases use the chosen result at each point (choices,
    from choices.toml; result A where none is chosen). base holds the tables
    as saved (all A); alternatives is the alternatives table.
    """
    label: str
    case: object
    summary: pd.DataFrame
    species: pd.DataFrame
    phases: pd.DataFrame
    comparison: pd.DataFrame | None = None
    path: Path | None = None
    alternatives: pd.DataFrame = field(default_factory=pd.DataFrame)
    choices: dict = field(default_factory=dict)
    base: tuple | None = None

    @classmethod
    def build(cls, label, case, summary, species, phases, comparison=None, path=None,
              alternatives=None, choices=None):
        """Tables with ``choices`` applied (the given tables are kept as ``base``)."""
        alternatives = pd.DataFrame() if alternatives is None else alternatives
        tables = cls(label, case, summary, species, phases, comparison, path, alternatives, {},
                     (summary, species, phases))
        return tables.with_choices(choices or {})

    def with_choices(self, choices):
        """These tables with another set of chosen results ({(T, P): label})."""
        base = self.base or (self.summary, self.species, self.phases)
        known = {(float(T), float(P)): label for (T, P), label in choices.items() if label != 'A'}
        summary, species, phases = apply_choices(*base, self.alternatives, known)
        return replace(self, summary=summary, species=species, phases=phases, choices=known, base=base)

    @property
    def flagged_points(self):
        """(T_C, P_bar) of the points with more than one result, in table order."""
        if self.alternatives.empty:
            return []
        points = self.alternatives[['T_C', 'P_bar']].drop_duplicates()
        return [(float(T), float(P)) for T, P in points.itertuples(index=False)]

    @property
    def pressures(self):
        return sorted(self.summary.P_bar.unique()) if 'P_bar' in self.summary else []

    @property
    def temperatures(self):
        return sorted(self.summary.T_C.unique()) if 'T_C' in self.summary else []

    @property
    def species_keys(self):
        """'phase:species' of every species with a nonzero amount."""
        if self.species.empty:
            return []
        part = self.species[self.species.mol_per_mol_C > 0]
        return sorted(set(part.phase + ':' + part.species))


def output_folder(results_dir, case_name, now=None):
    """results/<case name>-<date-time>, so runs never overwrite each other."""
    stamp = (now or datetime.now()).strftime('%Y%m%d-%H%M%S')
    return Path(results_dir) / f'{case_name}-{stamp}'


def saved_runs(results_dir):
    """Folders below results_dir that contain saved results, newest first."""
    root = Path(results_dir)
    if not root.exists():
        return []
    folders = [p.parent for p in root.glob('**/summary.csv')]
    return sorted(folders, key=lambda p: (p / 'summary.csv').stat().st_mtime, reverse=True)


def load_saved(folder, results_dir=None):
    folder = Path(folder)
    case = load_case(folder / 'case.toml') if (folder / 'case.toml').exists() else None
    comparison = None
    if (folder / 'comparison.csv').exists():
        comparison = pd.read_csv(folder / 'comparison.csv')
    label = str(folder.relative_to(results_dir)) if results_dir else folder.name
    summary, species, phases = (read_csv(folder / name) for name in ('summary.csv', 'species.csv', 'phases.csv'))
    if not phases.empty and 'P_bar' not in phases and 'P_bar' in summary:
        # results saved before pressure grids existed (one pressure)
        phases.insert(1, 'P_bar', summary.P_bar.iloc[0])
    alternatives = (read_csv(folder / 'alternatives.csv') if (folder / 'alternatives.csv').exists()
                    else pd.DataFrame())
    return RunTables.build(label, case, summary, species, phases, comparison, folder, alternatives,
                           read_choices(folder))


def live_tables(run):
    """Tables of a run in the pool (finished or still running)."""
    sweep = run.sweep_result()
    comparison = None
    if run.case.comparison and Path(run.case.comparison['file']).exists():
        comparison = pd.read_csv(run.case.comparison['file'])
    folder = getattr(run, 'output_dir', None)
    choices = read_choices(folder) if folder else {}
    return RunTables.build(f'{run.name} (this session, #{run.id})', run.case, sweep.summary(),
                           sweep.species_table(), sweep.phase_summary(), comparison, folder,
                           sweep.alternatives_table(), choices)
