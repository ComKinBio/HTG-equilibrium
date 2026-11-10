"""Plots of sweep results.

All functions work on the tables written by SweepResult.save (species.csv,
phases.csv, summary.csv), so saved results can be re-plotted without solving:

    htg plot results/my_case
"""
from pathlib import Path

import numpy as np
import pandas as pd

from .constants import PHASES

GROUPS = ['gas', 'aqueous neutrals', 'aqueous ions', 'organic liquid', 'solid']
QUANTITY_LABELS = {'mol_per_kg_dry_feed': 'mol / kg dry feed',
                   'mol_per_mol_C': 'mol / mol dry-feed C',
                   'g_per_kg_dry_feed': 'g / kg dry feed',
                   'carbon_mol_per_kg_dry_feed': 'mol C / kg dry feed',
                   'molality_mol_per_kg_water': 'mol / kg water'}


def _pyplot():
    import matplotlib.pyplot as plt
    return plt


def _grid(species, summary):
    if summary is not None and len(summary):
        return sorted(summary.T_C.astype(float).unique())
    return sorted(species.T_C.astype(float).unique())


def _solved(species, summary):
    if summary is not None and 'status' in summary:
        return set(summary.loc[summary.status == 'solved', 'T_C'].astype(float))
    return set(species.T_C.astype(float))


def species_profiles(species, quantity='mol_per_kg_dry_feed', summary=None, comparison=None,
                     comparison_label='Reference', top=12, min_value=1e-8, scale='symlog'):
    """One panel per phase group with the largest species against temperature.

    species     species table (SweepResult.species_table() or species.csv)
    comparison  wide table: T_C plus 'phase:species' columns, in units of ``quantity``
    top         species per panel (largest maximum first); None for all above min_value
    Failed temperatures appear as gaps.
    """
    plt = _pyplot()
    grid, solved = _grid(species, summary), _solved(species, summary)
    groups = [g for g in GROUPS if (species.group == g).any()]
    columns = 2 if len(groups) > 1 else 1
    rows = int(np.ceil(len(groups) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(6.5 * columns, 4 * rows), squeeze=False,
                             layout='constrained')
    placed = set()
    for ax, group in zip(axes.ravel(), groups):
        part = species[species.group == group]
        table = part.pivot_table(index='T_C', columns='species', values=quantity, aggfunc='sum')
        table = table.reindex(grid)
        for T in solved:
            if T in table.index:
                table.loc[T] = table.loc[T].fillna(0.0)
        peak = table.max().sort_values(ascending=False)
        names = [s for s in peak.index if peak[s] > min_value][:top]
        for name in names:
            line, = ax.plot(table.index, table[name], '-o', ms=3, label=name)
            key = f"{part.phase.iloc[0]}:{name}"
            if comparison is not None and key in comparison:
                ax.plot(comparison.T_C, comparison[key], 's', ms=4, mfc='none',
                        color=line.get_color())
                placed.add(key)
        ax.set(title=group, xlabel='Temperature [°C]', ylabel=QUANTITY_LABELS.get(quantity, quantity))
        if scale == 'symlog':
            ax.set_yscale('symlog', linthresh=max(min_value, 1e-6))
        ax.grid(alpha=0.25)
        if names:
            ax.legend(fontsize=7, ncols=2)
    for ax in axes.ravel()[len(groups):]:
        ax.axis('off')
    if comparison is not None:
        fig.suptitle(f'Lines: model; open squares: {comparison_label}')
        missing = [c for c in comparison.columns if c != 'T_C' and c not in placed]
        if missing:
            fig.text(0.01, 0.0, f'Not plotted (absent or below threshold in the model): {", ".join(missing)}',
                     fontsize=7, va='bottom')
    return fig


def carbon_and_mass(phases, summary=None):
    """Carbon allocation (% of feed carbon) and mass yield per phase against temperature."""
    plt = _pyplot()
    grid = sorted(summary.T_C.astype(float).unique()) if summary is not None and len(summary) \
        else sorted(phases.T_C.astype(float).unique())
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), layout='constrained')
    for ax, column, label in [(axes[0], 'carbon_percent', 'Carbon [% of feed C]'),
                              (axes[1], 'g_per_kg_dry_feed', 'Mass [g / kg dry feed]')]:
        table = phases.pivot_table(index='T_C', columns='phase', values=column).reindex(grid)
        for phase in [p for p in PHASES if p in table.columns]:
            if table[phase].abs().max() > 0:
                ax.plot(table.index, table[phase], '-o', ms=3, label=phase.replace('_', ' '))
        ax.set(xlabel='Temperature [°C]', ylabel=label)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    return fig


def read_csv(path):
    """A results CSV, or an empty table when the file is empty (no point solved)."""
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def save_plots(directory, quantity=None, comparison=None, comparison_label='Reference'):
    """Write species and carbon plots from the CSV files in ``directory``.

    With several pressures there is one pair of plots per pressure
    (species_240bar.png, carbon_240bar.png, ...).
    """
    directory = Path(directory)
    species = read_csv(directory / 'species.csv')
    phases = read_csv(directory / 'phases.csv')
    summary = read_csv(directory / 'summary.csv')
    if species.empty:
        return []
    pressures = sorted(summary.P_bar.unique()) if 'P_bar' in summary else [None]
    written = []
    for P in pressures:
        suffix = f'_{P:g}bar' if len(pressures) > 1 else ''
        if P is None:
            sp, ph, su = species, phases, summary
        else:
            sp, ph, su = species[species.P_bar == P], phases[phases.P_bar == P], summary[summary.P_bar == P]
        if sp.empty:
            continue
        fig = species_profiles(sp, quantity or 'mol_per_kg_dry_feed', su, comparison, comparison_label)
        if suffix:
            fig.suptitle(f'{P:g} bar' + (f' — {fig.get_suptitle()}' if fig.get_suptitle() else ''))
        fig.savefig(directory / f'species{suffix}.png', dpi=150)
        written.append(directory / f'species{suffix}.png')
        fig = carbon_and_mass(ph, su)
        if suffix:
            fig.suptitle(f'{P:g} bar')
        fig.savefig(directory / f'carbon{suffix}.png', dpi=150)
        written.append(directory / f'carbon{suffix}.png')
        _pyplot().close('all')
    return written
