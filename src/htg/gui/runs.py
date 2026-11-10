"""Result tables of runs, live (from the worker pool) or saved (results folders)."""
import tomllib
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import pandas as pd

from ..case import load_case
from ..plotting import read_csv


@dataclass
class RunTables:
    """Everything the result and comparison views need from one run."""
    label: str
    case: object
    summary: pd.DataFrame
    species: pd.DataFrame
    phases: pd.DataFrame
    comparison: pd.DataFrame | None = None
    path: Path | None = None

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
    return RunTables(label, case, summary, species, phases, comparison, folder)


def live_tables(run):
    """Tables of a run in the pool (finished or still running)."""
    sweep = run.sweep_result()
    comparison = None
    if run.case.comparison and Path(run.case.comparison['file']).exists():
        comparison = pd.read_csv(run.case.comparison['file'])
    return RunTables(f'{run.name} (this session, #{run.id})', run.case, sweep.summary(),
                     sweep.species_table(), sweep.phase_summary(), comparison, run.output_dir)


def flatten(data, prefix=''):
    """{'feed': {'water_g_per_g_dry': 4}} -> {'feed.water_g_per_g_dry': 4}."""
    out = {}
    for key, value in data.items():
        name = f'{prefix}{key}'
        if isinstance(value, dict) and value:
            out.update(flatten(value, name + '.'))
        else:
            out[name] = value
    return out


def case_differences(cases):
    """Settings that differ between cases: DataFrame with one column per label.

    cases  {label: Case}
    """
    flat = {}
    for label, case in cases.items():
        if case is None:
            continue
        data = tomllib.loads(case.to_toml())
        data.pop('name', None)
        data.pop('description', None)
        data.get('comparison', {}).pop('file', None)
        flat[label] = flatten(data)
    keys = sorted(set().union(*[set(f) for f in flat.values()])) if flat else []
    rows = []
    for key in keys:
        values = [flat[label].get(key) for label in flat]
        if any(v != values[0] for v in values[1:]):
            rows.append({'setting': key, **{label: _show(v) for label, v in zip(flat, values)}})
    return pd.DataFrame(rows, columns=['setting', *flat])


def _show(value):
    if value is None:
        return '(default)'
    if isinstance(value, list):
        return ', '.join(map(str, value)) if len(value) <= 12 else f'{len(value)} values'
    return str(value)
