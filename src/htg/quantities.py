"""Quantities shown on the Results tab, computed from the saved result tables.

The tables written by a run (species.csv, phases.csv, summary.csv) give amounts
per kg of dry feed (solids including ash). With

    s = solids / total feed      (total feed = solids + water)
    a = ash / solids

(htg.feed_v2.feed_fractions gives s and a of a run's feed), the other bases
follow by plain scaling:

    per kg total feed          = per kg dry feed * s
    per kg dry, ash-free feed  = per kg dry feed / (1 - a)

Every function here is pure (pandas in, pandas out) and needs neither the GUI
nor Julia. Long tables returned by the amount functions have the columns
T_C, P_bar, [phase,] species, value.
"""
import numpy as np
import pandas as pd

WATER = 'water'
TINY = 1e-8  # amounts at or below this (mol / mol C) count as "none" (solver noise is ~1e-10)

BASES = {'total': 'total feed', 'dry': 'dry feed', 'daf': 'dry ash-free feed'}
BASIS_OPTIONS = {'total': 'total feed (solids + water)', 'dry': 'dry feed (solids incl. ash)',
                 'daf': 'dry, ash-free feed'}
AMOUNTS = {'mol': 'mol', 'g': 'g'}
NO_TOTAL_MASS = 'This run has no total mass.'

_AMOUNT_COLUMNS = {'mol': 'mol_per_kg_dry_feed', 'g': 'g_per_kg_dry_feed'}


def basis_factor(basis, s=None, a=0.0):
    """Factor from 'per kg dry feed' to the basis 'total', 'dry' or 'daf'."""
    if basis == 'dry':
        return 1.0
    if basis == 'total':
        if s is None:
            raise ValueError(NO_TOTAL_MASS)
        return float(s)
    if basis == 'daf':
        a = float(a or 0.0)
        if not 0 <= a < 1:
            raise ValueError(f'ash fraction must be in [0, 1), not {a}')
        return 1.0 / (1.0 - a)
    raise ValueError(f'unknown basis {basis!r}')


def available_bases(s):
    """{basis: None if it can be used, else the tooltip saying why not}."""
    return {'total': None if s is not None else NO_TOTAL_MASS, 'dry': None, 'daf': None}


def amount_unit(amount='mol', basis='dry'):
    """'mol / kg total feed' and similar."""
    return f'{AMOUNTS[amount]} / kg {BASES[basis]}'


# --- amounts ------------------------------------------------------------------------------------

def amounts(species, amount='mol', basis='dry', s=None, a=0.0):
    """Amount of every (phase, species) in the chosen unit.

    Long table T_C, P_bar, phase, species, value.
    """
    factor = basis_factor(basis, s, a)
    out = species[['T_C', 'P_bar', 'phase', 'species']].copy()
    out['value'] = species[_AMOUNT_COLUMNS[amount]].astype(float) * factor
    return out.reset_index(drop=True)


def total_amounts(species, amount='mol', basis='dry', s=None, a=0.0):
    """Amount of each species summed over all phases: T_C, P_bar, species, value."""
    part = amounts(species, amount, basis, s, a)
    return part.groupby(['T_C', 'P_bar', 'species'], as_index=False, sort=False)['value'].sum()


def gas_composition(species, by='mol', exclude_water=False):
    """Gas-phase composition in mol % or mass %: T_C, P_bar, species, value.

    Normalised over the gas phase at each point; with exclude_water the water is
    removed before normalising (dry gas).
    """
    gas = species[species.phase == 'gas']
    if exclude_water:
        gas = gas[gas.species != WATER]
    column = _AMOUNT_COLUMNS['mol' if by == 'mol' else 'g']
    out = gas[['T_C', 'P_bar', 'species']].copy()
    values = gas[column].astype(float)
    totals = values.groupby([gas.T_C, gas.P_bar]).transform('sum')
    out['value'] = np.where(totals > 0, 100.0 * values / totals.where(totals > 0, 1.0), np.nan)
    return out.reset_index(drop=True)


def carbon_gasification_efficiency(summary):
    """Carbon in the gas phase as % of feed carbon (the saved gas_carbon_percent).

    Only gas-phase carbon counts. T_C, P_bar, value; only solved points.
    """
    part = solved_rows(summary)
    out = part[['T_C', 'P_bar']].copy()
    out['value'] = part['gas_carbon_percent'].astype(float)
    return out.reset_index(drop=True)


def phase_carbon(phases, summary=None):
    """Share of feed carbon in each phase [%]: T_C, P_bar, phase, value."""
    phases = with_pressure(phases, summary)
    out = phases[['T_C', 'P_bar', 'phase']].copy()
    out['value'] = phases['carbon_percent'].astype(float)
    return out.reset_index(drop=True)


def with_pressure(table, summary):
    """Add P_bar to tables saved before pressure grids existed (one pressure)."""
    if 'P_bar' in table or summary is None or 'P_bar' not in summary or summary.empty:
        return table
    table = table.copy()
    table.insert(1, 'P_bar', float(summary.P_bar.iloc[0]))
    return table


# --- profiles along one axis --------------------------------------------------------------------

def solved_rows(summary):
    return summary[summary.status == 'solved'] if 'status' in summary else summary


def at(table, column, value):
    """Rows where column is (numerically) equal to value."""
    return table[np.isclose(table[column].astype(float), float(value), rtol=1e-9, atol=1e-9)]


def grid_values(summary, axis):
    """Sorted values of 'T_C' or 'P_bar' in the run."""
    return sorted(float(v) for v in summary[axis].unique()) if axis in summary else []


def profile(long, summary, axis='T_C', fixed=None, keys=None, key_column='species'):
    """Values along one axis at a fixed value of the other.

    long        table with T_C, P_bar, key_column and value
    axis        'T_C' or 'P_bar'; the other one is held at ``fixed``
    keys        columns to return (default: every key in ``long``)
    Returns a DataFrame indexed by the axis values of every point in the summary at
    ``fixed``. Failed (or unfinished) points are NaN; a key missing at a solved point
    is 0.
    """
    other = 'P_bar' if axis == 'T_C' else 'T_C'
    rows = at(summary, other, fixed) if fixed is not None else summary
    grid = sorted(float(v) for v in rows[axis].unique())
    solved = {float(v) for v in solved_rows(rows)[axis]}
    part = at(long, other, fixed) if fixed is not None else long
    if keys is None:
        keys = list(dict.fromkeys(part[key_column]))
    if part.empty:
        table = pd.DataFrame(index=pd.Index(grid, dtype=float), columns=list(keys), dtype=float)
    else:
        table = part.pivot_table(index=axis, columns=key_column, values='value', aggfunc='sum')
        table.index = table.index.astype(float)
        table = table.reindex(index=grid, columns=list(keys)).astype(float)
    for value in grid:
        if value in solved:
            table.loc[value] = table.loc[value].fillna(0.0)
        else:
            table.loc[value] = np.nan
    table.index.name = axis
    return table


def largest(long, n=8, include_water=False, key_column='species'):
    """The n keys with the largest maximum value (water left out unless asked for)."""
    part = long if include_water else long[long[key_column] != WATER]
    if part.empty:
        return []
    peak = part.groupby(key_column)['value'].max().sort_values(ascending=False, kind='stable')
    return [k for k in peak.index if peak[k] > 0][:n]


def phase_is_empty(species, phase, points=None, threshold=TINY):
    """True when the phase has no amount (mol / mol C > threshold) at any of the points.

    points  optional DataFrame/iterable of (T_C, P_bar) to look at; default all points
    """
    part = species[species.phase == phase]
    if points is not None:
        keys = {(float(T), float(P)) for T, P in points}
        mask = [(float(T), float(P)) in keys for T, P in zip(part.T_C, part.P_bar)]
        part = part[np.asarray(mask, dtype=bool)] if len(part) else part
    return part.empty or not (part.mol_per_mol_C.astype(float) > threshold).any()


# --- measured data ------------------------------------------------------------------------------

def measured_amounts(table, quantity, species, amount='mol', basis='dry', s=None, a=0.0):
    """Measured (reference) data converted to the plotted amount unit.

    table     wide table: T_C (optionally P_bar) plus 'phase:species' columns, as in the
              comparison.csv files of a case
    quantity  unit of the table: 'mol_per_mol_C', 'mol_per_kg_dry_feed' or 'g_per_kg_dry_feed'
    species   the run's species table (gives the feed-carbon factor and molar masses)
    Long table T_C, [P_bar,] phase, species, value; values that cannot be converted are
    left out.
    """
    factor = basis_factor(basis, s, a)
    per_C = feed_carbon_per_kg(species)
    masses = molar_masses(species)
    id_columns = [c for c in ('T_C', 'P_bar') if c in table]
    rows = []
    for column in table.columns:
        if column in id_columns or ':' not in str(column):
            continue
        phase, name = str(column).split(':', 1)
        values = table[column].astype(float)
        if quantity == 'mol_per_mol_C':
            if per_C is None:
                continue
            values = values * per_C
        elif quantity == 'g_per_kg_dry_feed' and amount == 'mol':
            if name not in masses:
                continue
            values = values / masses[name]
        elif quantity not in ('mol_per_kg_dry_feed', 'g_per_kg_dry_feed'):
            continue
        if amount == 'g' and quantity != 'g_per_kg_dry_feed':
            if name not in masses:
                continue
            values = values * masses[name]
        part = table[id_columns].copy()
        part['phase'], part['species'], part['value'] = phase, name, values * factor
        rows.append(part.dropna(subset=['value']))
    columns = [*id_columns, 'phase', 'species', 'value']
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=columns)


def feed_carbon_per_kg(species):
    """mol of feed carbon per kg dry feed (= mol_per_kg_dry_feed / mol_per_mol_C), or None."""
    part = species[species.mol_per_mol_C.astype(float) > TINY]
    if part.empty:
        return None
    return float(np.median(part.mol_per_kg_dry_feed.astype(float) / part.mol_per_mol_C.astype(float)))


def molar_masses(species):
    """{species name: g/mol} from the species table."""
    part = species[species.mol_per_kg_dry_feed.astype(float) > 0]
    ratio = part.g_per_kg_dry_feed.astype(float) / part.mol_per_kg_dry_feed.astype(float)
    return {name: float(v) for name, v in ratio.groupby(part.species).median().items()}
