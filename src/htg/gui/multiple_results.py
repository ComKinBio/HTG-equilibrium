"""Points with more than one result: banner texts, the dialog and comparison mode.

Different starting points can end at different minima (results A, B, C, … by
increasing Gibbs energy; htg.results.Alternative). The alternatives table of a run
(alternatives.csv, RunTables.alternatives) lists every result at such points. The
dialog lets the user choose the result each point uses in plots and tables
(written to choices.toml) and which results are drawn together in comparison mode.
"""
import numpy as np
import pandas as pd
from nicegui import ui

from .. import units
from ..plotbuild import PHASE_TITLES
from ..plotspec import RESULT_DASHES

DASH_WORDS = {'solid': 'solid', 'dash': 'dashed', 'dot': 'dotted', 'dashdot': 'dash-dot'}
CAPTION = ('Different starting points ended at different minima. Plots and tables use the lowest Gibbs '
           'energy unless you choose another result. Compare in plots draws the ticked results together, '
           'one line style each.')
TITLE = 'Points with more than one result'
COLUMNS = ('Result', 'Gibbs energy [kJ / mol C]', 'Above lowest [J / mol C]', 'Gas carbon [%]',
           'Starts that found it', 'Use in plots and tables', 'Show when comparing')
LARGEST_N = 5


# ---------------------------------------------------------------- pure helpers
def _at(table, point):
    T, P = point
    mask = np.isclose(table.T_C.astype(float), float(T))
    if 'P_bar' in table:
        mask &= np.isclose(table.P_bar.astype(float), float(P))
    return table[mask]


def _join(words):
    return words[0] if len(words) == 1 else ', '.join(words[:-1]) + ' and ' + words[-1]


def banner_text(n_flagged, n_points):
    """'3 of 13 points ended at more than one result (…).'"""
    return (f'{n_flagged} of {n_points} point{"s" if n_points != 1 else ""} ended at more than one result '
            '(different starting points gave different minima).')


def comparing_text(n_points, labels):
    """'Comparing results at 3 points: A (solid, lowest Gibbs energy) and B (dashed). Elsewhere …'."""
    parts = []
    for i, label in enumerate(labels):
        dash = DASH_WORDS[RESULT_DASHES[i % len(RESULT_DASHES)]]
        parts.append(f'{label} ({dash}, lowest Gibbs energy)' if i == 0 else f'{label} ({dash})')
    return (f'Comparing results at {n_points} point{"s" if n_points != 1 else ""}: {_join(parts)}. '
            'Elsewhere there is only one result.')


def point_text(point, t_unit='C', p_unit='bar'):
    """'400 °C · 240 bar' in the display units."""
    T, P = point
    return f'{units.display(T, t_unit)} {units.label(t_unit)} · {units.display(P, p_unit)} {units.label(p_unit)}'


def point_results(alternatives, point):
    """One row per result at a point: result, gibbs_kJ, above_J, gas_carbon, n_starts, starts."""
    rows = _at(alternatives, point)
    out = []
    if rows.empty:
        return pd.DataFrame(columns=['result', 'gibbs_kJ', 'above_J', 'gas_carbon', 'n_starts', 'starts'])
    lowest = float(rows.gibbs_J_per_mol_C.min())
    for label, part in rows.groupby('result', sort=True):
        gibbs = float(part.gibbs_J_per_mol_C.iloc[0])
        gas_carbon = np.nan
        if 'carbon_mol_per_kg_dry_feed' in part:
            carbon = part.carbon_mol_per_kg_dry_feed.astype(float)
            total = float(carbon.sum())
            gas_carbon = 100.0 * float(carbon[part.phase == 'gas'].sum()) / total if total > 0 else 0.0
        starts = str(part.start_labels.iloc[0]) if 'start_labels' in part and pd.notna(part.start_labels.iloc[0]) \
            else ''
        out.append({'result': str(label), 'gibbs_kJ': gibbs / 1000.0, 'above_J': gibbs - lowest,
                    'gas_carbon': gas_carbon, 'n_starts': int(part.n_starts.iloc[0]), 'starts': starts})
    return pd.DataFrame(out)


def largest_differences(alternatives, point, factor=1.0, n=LARGEST_N):
    """The n (phase, species) whose amounts differ most between the results at a point.

    Amounts are mol / kg dry feed times ``factor`` (s for per kg total feed). Columns:
    'Species' ('Gas · methane'), then one column per result.
    """
    rows = _at(alternatives, point)
    if rows.empty:
        return pd.DataFrame(columns=['Species'])
    column = 'mol_per_kg_dry_feed' if 'mol_per_kg_dry_feed' in rows else 'mol_per_mol_C'
    table = rows.pivot_table(index=['phase', 'species'], columns='result', values=column, aggfunc='sum')
    table = table.fillna(0.0) * float(factor)
    spread = table.max(axis=1) - table.min(axis=1)
    order = spread.sort_values(ascending=False, kind='stable').index[:n]
    out = table.loc[order].reset_index()
    out.insert(0, 'Species', [f'{PHASE_TITLES.get(p, p)} · {s}' for p, s in zip(out.phase, out.species)])
    return out.drop(columns=['phase', 'species']).rename_axis(columns=None)


def results_at(alternatives, point):
    rows = _at(alternatives, point)
    return sorted(str(r) for r in rows.result.unique()) if len(rows) else []


# ---------------------------------------------------------------- the dialog
class MultipleResultsDialog:
    """The "Points with more than one result" dialog.

    tables       RunTables of the run (alternatives, choices, flagged_points)
    show         {(T, P): {label: bool}}: which results comparison mode draws (changed in place)
    on_apply(choices)   called by "Apply choices" with {(T, P): label}
    on_compare()        called by "Compare in plots"
    factor, unit_text   amounts in the "Largest differences" list: mol / kg dry feed × factor
    """

    def __init__(self, tables, show, on_apply, on_compare, t_unit='C', p_unit='bar', factor=1.0,
                 unit_text='mol / kg total feed'):
        self.tables = tables
        self.alternatives = tables.alternatives
        self.points = tables.flagged_points
        self.choices = {point: tables.choices.get(point, 'A') for point in self.points}
        self.show = show
        for point in self.points:
            entry = show.setdefault(point, {})
            for label in results_at(self.alternatives, point):
                entry.setdefault(label, True)
        self.on_apply, self.on_compare = on_apply, on_compare
        self.t_unit, self.p_unit = t_unit, p_unit
        self.factor, self.unit_text = factor, unit_text
        self.current = self.points[0] if self.points else None
        self.items = {}
        with ui.dialog() as self.dialog, ui.card().classes('w-[64rem] max-w-full'):
            ui.label(TITLE).classes('text-h6')
            ui.label(CAPTION).classes('text-body2 text-grey-8')
            with ui.row().classes('w-full no-wrap gap-4'):
                with ui.list().props('bordered separator dense').classes('w-56 shrink-0 max-h-96 overflow-auto'):
                    for i, point in enumerate(self.points):
                        item = ui.item(point_text(point, t_unit, p_unit), on_click=lambda p=point: self.pick(p)) \
                            .mark(f'point-{i}')
                        self.items[point] = item
                self.detail = ui.column().classes('grow gap-2')
            with ui.row().classes('w-full justify-end gap-2'):
                ui.button('Close', on_click=self.dialog.close).props('flat no-caps')
                ui.button('Compare in plots', on_click=self.compare).props('outline no-caps')
                ui.button('Apply choices', on_click=self.apply).props('no-caps')
        if self.current is not None:
            self.pick(self.current)

    def pick(self, point):
        self.current = point
        for p, item in self.items.items():
            item.classes(remove='bg-blue-1 text-weight-medium')
            if p == point:
                item.classes(add='bg-blue-1 text-weight-medium')
        self.detail.clear()
        results = point_results(self.alternatives, point)
        self.radios, self.checks = {}, {}
        with self.detail:
            with ui.grid(columns=len(COLUMNS)).classes('w-full items-center gap-x-4 gap-y-1'):
                for name in COLUMNS:
                    ui.label(name).classes('text-caption text-weight-medium')
                for row in results.itertuples(index=False):
                    label = row.result
                    ui.label(label)
                    ui.label(units.format_number(row.gibbs_kJ))
                    ui.label(units.format_number(row.above_J, digits=4))
                    ui.label('' if pd.isna(row.gas_carbon) else units.format_number(row.gas_carbon, digits=4))
                    starts = ui.label(str(row.n_starts))
                    if row.starts:
                        starts.tooltip(row.starts)
                    self.radios[label] = ui.radio({label: ''}, value=label if self.choices[point] == label else None,
                                                  on_change=lambda e, lb=label: self.choose(lb, e.value)) \
                        .props('dense')
                    check = ui.checkbox(value=self.show[point].get(label, True),
                                        on_change=lambda e, lb=label: self.tick(lb, e.value)).props('dense')
                    if label == 'A':  # result A is always drawn
                        check.set_value(True)
                        check.disable()
                    self.checks[label] = check
            ui.label(f'Largest differences [{self.unit_text}]').classes('text-weight-medium q-mt-sm')
            differences = largest_differences(self.alternatives, point, self.factor)
            with ui.grid(columns=len(differences.columns)).classes('gap-x-4 gap-y-0'):
                for name in differences.columns:
                    ui.label(name).classes('text-caption text-weight-medium')
                for row in differences.itertuples(index=False):
                    for i, value in enumerate(row):
                        ui.label(value if i == 0 else units.format_number(value, digits=4))

    def choose(self, label, value):
        if value != label:
            return
        self.choices[self.current] = label
        for other, radio in self.radios.items():
            if other != label and radio.value is not None:
                radio.set_value(None)

    def tick(self, label, value):
        self.show[self.current][label] = bool(value)

    def apply(self):
        self.dialog.close()
        self.on_apply(dict(self.choices))

    def compare(self):
        self.dialog.close()
        self.on_compare()

    def open(self):
        self.dialog.open()
        return self
