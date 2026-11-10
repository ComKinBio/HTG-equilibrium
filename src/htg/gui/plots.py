"""Interactive Plotly figures of result tables.

Every function takes the tables of htg.results (species, phases, summary) and
returns a plotly Figure. Failed points appear as gaps; species absent at a
solved point count as zero.
"""
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from ..constants import PHASES

GROUPS = ['gas', 'aqueous neutrals', 'aqueous ions', 'organic liquid', 'solid']
QUANTITIES = {
    'mol_per_kg_dry_feed': 'mol / kg dry feed',
    'mol_per_mol_C': 'mol / mol dry-feed C',
    'g_per_kg_dry_feed': 'g / kg dry feed',
    'carbon_mol_per_kg_dry_feed': 'mol C / kg dry feed',
    'molality_mol_per_kg_water': 'mol / kg water (ions)',
}
COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
          '#bcbd22', '#17becf', '#393b79', '#ad494a', '#637939', '#7b4173', '#3182bd', '#e6550d']
DASHES = ['solid', 'dash', 'dot', 'dashdot', 'longdash', 'longdashdot']
PHASE_LABELS = {'gas': 'gas', 'aqueous': 'aqueous', 'organic_liquid': 'organic liquid', 'solid': 'solid'}


def _layout(fig, title, x_title, y_title, log=False):
    fig.update_layout(title=title, xaxis_title=x_title, yaxis_title=y_title, template='plotly_white',
                      legend={'font': {'size': 11}}, margin={'l': 60, 'r': 20, 't': 50, 'b': 50},
                      hovermode='closest')
    if log:
        fig.update_yaxes(type='log', exponentformat='power')
    return fig


def _solved(summary, P):
    part = summary[(summary.P_bar == P)] if 'P_bar' in summary else summary
    if 'status' in part:
        part = part[part.status == 'solved']
    return set(part.T_C)


def _profile(species, summary, P, key_column, keys, quantity, axis='T_C', fixed=None):
    """Table index = axis values, columns = keys, gaps at failed points."""
    if axis == 'T_C':
        part = species[species.P_bar == P]
        grid = sorted(summary.loc[summary.P_bar == P, 'T_C'].unique())
        solved = _solved(summary, P)
    else:
        part = species[species.T_C == fixed]
        rows = summary[summary.T_C == fixed]
        grid = sorted(rows.P_bar.unique())
        solved = set(rows.loc[rows.status == 'solved', 'P_bar']) if 'status' in rows else set(grid)
    table = part.pivot_table(index=axis, columns=key_column, values=quantity, aggfunc='sum')
    table = table.reindex(index=grid, columns=keys)
    for value in solved:
        if value in table.index:
            table.loc[value] = table.loc[value].fillna(0.0)
    return table


def _top(species, group, quantity, top, min_value):
    part = species[species.group == group]
    if part.empty:
        return []
    peak = part.groupby('species')[quantity].max().sort_values(ascending=False)
    names = [s for s in peak.index if peak[s] > min_value]
    return names[:top] if top else names


def _for_log(values, log, floor):
    values = np.asarray(values, dtype=float)
    return np.where(values > floor, values, np.nan) if log else values


def species_vs_temperature(species, summary, group, quantity='mol_per_kg_dry_feed', pressures=None,
                           top=12, log=True, min_value=1e-8, comparison=None, comparison_label='Reference'):
    """Species of one group against temperature; colour = species, line style = pressure."""
    fig = go.Figure()
    if species.empty:
        return _layout(fig, group, 'Temperature [°C]', QUANTITIES.get(quantity, quantity), log)
    pressures = pressures or sorted(summary.P_bar.unique())
    names = _top(species[species.P_bar.isin(pressures)], group, quantity, top, min_value)
    part = species[species.group == group]
    phase = part.phase.iloc[0] if len(part) else None
    for k, P in enumerate(pressures):
        table = _profile(part, summary, P, 'species', names, quantity)
        for i, name in enumerate(names):
            label = name if len(pressures) == 1 else f'{name}, {P:g} bar'
            fig.add_trace(go.Scatter(
                x=table.index, y=_for_log(table[name], log, min_value), name=label, mode='lines+markers',
                legendgroup=name, line={'color': COLORS[i % len(COLORS)], 'dash': DASHES[k % len(DASHES)]},
                marker={'size': 5},
                hovertemplate=f'{name}<br>%{{x:g}} °C, {P:g} bar<br>%{{y:.4g}}<extra></extra>'))
    if comparison is not None and phase is not None:
        for i, name in enumerate(names):
            key = f'{phase}:{name}'
            if key in comparison:
                fig.add_trace(go.Scatter(
                    x=comparison.T_C, y=_for_log(comparison[key], log, min_value), mode='markers',
                    name=f'{name} ({comparison_label})', legendgroup=name,
                    marker={'symbol': 'square-open', 'size': 8, 'color': COLORS[i % len(COLORS)]}))
    return _layout(fig, group, 'Temperature [°C]', QUANTITIES.get(quantity, quantity), log)


def species_vs_pressure(species, summary, group, temperature, quantity='mol_per_kg_dry_feed', top=12,
                        log=True, min_value=1e-8):
    """Species of one group against pressure at one temperature."""
    fig = go.Figure()
    part = species[(species.group == group) & (species.T_C == temperature)]
    names = _top(part, group, quantity, top, min_value)
    if names:
        table = _profile(part, summary, None, 'species', names, quantity, axis='P_bar', fixed=temperature)
        for i, name in enumerate(names):
            fig.add_trace(go.Scatter(x=table.index, y=_for_log(table[name], log, min_value), name=name,
                                     mode='lines+markers', line={'color': COLORS[i % len(COLORS)]}))
    return _layout(fig, f'{group} at {temperature:g} °C', 'Pressure [bar]',
                   QUANTITIES.get(quantity, quantity), log)


def map_items(tables):
    """Choices for the T-P map: summary quantities, phase carbon shares, species."""
    items = {'gas_carbon_percent': 'Gas carbon [% of feed C]',
             'gibbs_J_per_mol_C': 'Gibbs energy [J / mol C]'}
    items.update({f'carbon:{p}': f'Carbon in {PHASE_LABELS[p]} [% of feed C]' for p in PHASES})
    items.update({f'species:{key}': key.replace(':', ': ') for key in tables.species_keys})
    return items


def tp_map(tables, item, quantity='mol_per_kg_dry_feed'):
    """Colour map of one quantity over the temperature-pressure grid."""
    summary = tables.summary
    if item.startswith('species:'):
        phase, name = item[len('species:'):].split(':', 1)
        part = tables.species[(tables.species.phase == phase) & (tables.species.species == name)]
        values = part.pivot_table(index='P_bar', columns='T_C', values=quantity, aggfunc='sum')
        title, unit = f'{name} ({PHASE_LABELS.get(phase, phase)})', QUANTITIES.get(quantity, quantity)
        fill = 0.0
    elif item.startswith('carbon:'):
        phase = item[len('carbon:'):]
        part = tables.phases[tables.phases.phase == phase]
        values = part.pivot_table(index='P_bar', columns='T_C', values='carbon_percent', aggfunc='sum')
        title, unit, fill = f'Carbon in {PHASE_LABELS[phase]}', '% of feed C', 0.0
    else:
        part = summary[summary.status == 'solved'] if 'status' in summary else summary
        values = part.pivot_table(index='P_bar', columns='T_C', values=item, aggfunc='first')
        title, unit, fill = map_items(tables).get(item, item), '', None
    values = values.reindex(index=tables.pressures, columns=tables.temperatures)
    if fill is not None:  # solved points without the species are zero, failed points stay empty
        for _, row in summary[summary.status == 'solved'].iterrows():
            if pd.isna(values.loc[row.P_bar, row.T_C]):
                values.loc[row.P_bar, row.T_C] = fill
    fig = go.Figure(go.Heatmap(z=values.values, x=[f'{t:g}' for t in values.columns],
                               y=[f'{p:g}' for p in values.index], colorscale='Viridis',
                               colorbar={'title': unit},
                               hovertemplate='%{x} °C, %{y} bar<br>%{z:.4g}<extra></extra>'))
    fig.update_xaxes(type='category')
    fig.update_yaxes(type='category')
    return _layout(fig, title, 'Temperature [°C]', 'Pressure [bar]')


def carbon_vs_temperature(phases, summary, pressures=None):
    """Share of feed carbon in each phase against temperature."""
    fig = go.Figure()
    if phases.empty:
        return _layout(fig, 'Carbon allocation', 'Temperature [°C]', '% of feed carbon')
    pressures = pressures or sorted(summary.P_bar.unique())
    for k, P in enumerate(pressures):
        table = _profile(phases, summary, P, 'phase', list(PHASES), 'carbon_percent')
        for i, phase in enumerate(PHASES):
            if table[phase].abs().max() > 0:
                label = PHASE_LABELS[phase] if len(pressures) == 1 else f'{PHASE_LABELS[phase]}, {P:g} bar'
                fig.add_trace(go.Scatter(x=table.index, y=table[phase], name=label, mode='lines+markers',
                                         legendgroup=phase,
                                         line={'color': COLORS[i], 'dash': DASHES[k % len(DASHES)]}))
    return _layout(fig, 'Carbon allocation', 'Temperature [°C]', '% of feed carbon')


def compare_runs(runs, key, quantity='mol_per_kg_dry_feed', log=False, comparison=None,
                 comparison_label='Reference', min_value=1e-8):
    """One species ('phase:species') or phase carbon share ('carbon:phase') across runs."""
    fig = go.Figure()
    for i, tables in enumerate(runs):
        for k, P in enumerate(tables.pressures):
            if key.startswith('carbon:'):
                phase = key[len('carbon:'):]
                table = _profile(tables.phases, tables.summary, P, 'phase', [phase], 'carbon_percent')
                y = table[phase]
            else:
                phase, name = key.split(':', 1)
                part = tables.species[tables.species.phase == phase]
                table = _profile(part, tables.summary, P, 'species', [name], quantity)
                y = table[name]
            label = tables.label if len(tables.pressures) == 1 else f'{tables.label}, {P:g} bar'
            fig.add_trace(go.Scatter(x=table.index, y=_for_log(y, log, min_value), name=label,
                                     mode='lines+markers',
                                     line={'color': COLORS[i % len(COLORS)], 'dash': DASHES[k % len(DASHES)]}))
    if comparison is not None and key in comparison:
        fig.add_trace(go.Scatter(x=comparison.T_C, y=_for_log(comparison[key], log, min_value),
                                 name=comparison_label, mode='markers',
                                 marker={'symbol': 'square-open', 'size': 9, 'color': 'black'}))
    unit = '% of feed carbon' if key.startswith('carbon:') else QUANTITIES.get(quantity, quantity)
    return _layout(fig, key.replace(':', ': '), 'Temperature [°C]', unit, log)
