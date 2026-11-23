"""Points with more than one result (WP6): banner, dialog, choices.toml and comparison mode,
with NiceGUI's simulated user on a stand-in run (no Julia)."""
import asyncio
import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from nicegui import ui
from test_plot_core import AMOUNTS, make_tables

from htg.gui import multiple_results as mr
from htg.gui.results_tab import ResultsTab
from htg.gui.runs import load_saved
from htg.results import read_choices
from htg.settings import Settings

TEMPLATE = Path(__file__).resolve().parents[1] / 'src/htg/templates/case_template.toml'


def _result(species, T, label, gibbs, n_starts, changes=None):
    """Alternative rows at T: the species rows there, with amounts scaled by ``changes``."""
    rows = species[species.T_C == T].copy()
    for (phase, name), factor in (changes or {}).items():
        mask = (rows.phase == phase) & (rows.species == name)
        for column in ('mol_per_mol_C', 'mol_per_kg_dry_feed', 'g_per_kg_dry_feed', 'carbon_mol_per_kg_dry_feed'):
            rows.loc[mask, column] = rows.loc[mask, column] * factor
    rows.insert(2, 'result', label)
    rows.insert(3, 'gibbs_J_per_mol_C', gibbs)
    rows.insert(4, 'n_starts', n_starts)
    rows['phase_gibbs_J_per_mol_C'] = np.nan
    rows['start_labels'] = f'start {label}'
    return rows


def stand_in_run(folder):
    """make_tables (240 bar; 300, 400, 450 failed, 500, 600, 700 °C) saved with alternatives at 400 (A, B)
    and 500 °C (A, B, C)."""
    tables = make_tables({**AMOUNTS, 600: AMOUNTS[500], 700: AMOUNTS[500]})
    species = tables.species
    alternatives = pd.concat([
        _result(species, 400.0, 'A', -1e6, 3),
        _result(species, 400.0, 'B', -1e6 + 120.0, 2, {('gas', 'carbon dioxide'): 1.5, ('solid', 'carbon'): 0.5}),
        _result(species, 500.0, 'A', -1e6, 2),
        _result(species, 500.0, 'B', -1e6 + 30.0, 1, {('gas', 'methane'): 2.0}),
        _result(species, 500.0, 'C', -1e6 + 900.0, 1, {('gas', 'hydrogen'): 3.0}),
    ], ignore_index=True)
    first = ['T_C', 'P_bar', 'result', 'gibbs_J_per_mol_C', 'n_starts', 'phase', 'species', 'mol_per_mol_C']
    alternatives = alternatives[first + [c for c in alternatives.columns if c not in first]]
    folder.mkdir(parents=True, exist_ok=True)
    tables.summary.to_csv(folder / 'summary.csv', index=False)
    species.to_csv(folder / 'species.csv', index=False)
    tables.phases.to_csv(folder / 'phases.csv', index=False)
    alternatives.to_csv(folder / 'alternatives.csv', index=False)
    shutil.copyfile(TEMPLATE, folder / 'case.toml')  # format 2, with a total mass
    return folder


class SavedState:
    """Reads the run from disk every time (as htg.gui.app.State does for saved runs)."""

    def __init__(self, folder, settings):
        self.folder, self.settings = folder, settings

    def sources(self):
        return {'saved:run': 'run'}

    def run(self, source):
        return None

    def tables(self, source):
        return load_saved(self.folder, self.folder.parent) if source == 'saved:run' else None


def test_texts_and_tables(tmp_path):
    tables = load_saved(stand_in_run(tmp_path / 'run'), tmp_path)
    assert tables.flagged_points == [(400.0, 240.0), (500.0, 240.0)]
    assert mr.banner_text(3, 13) == ('3 of 13 points ended at more than one result (different starting points '
                                     'gave different minima).')
    assert mr.comparing_text(3, ['A', 'B']) == ('Comparing results at 3 points: A (solid, lowest Gibbs energy) '
                                                'and B (dashed). Elsewhere there is only one result.')
    assert mr.comparing_text(1, ['A', 'B', 'C']).startswith(
        'Comparing results at 1 point: A (solid, lowest Gibbs energy), B (dashed) and C (dotted).')
    assert mr.point_text((400.0, 240.0), 'K', 'MPa') == '673.15 K · 24 MPa'
    rows = mr.point_results(tables.alternatives, (400.0, 240.0)).set_index('result')
    assert list(rows.index) == ['A', 'B']
    assert rows.loc['B', 'gibbs_kJ'] == pytest.approx(-999.88) and rows.loc['B', 'above_J'] == pytest.approx(120.0)
    assert rows.loc['A', 'gas_carbon'] == pytest.approx(100 * 30 / 40)  # CO2 20 + CH4 10 of 40 mol C
    assert rows.loc['B', 'gas_carbon'] == pytest.approx(100 * 40 / 46)  # CO2 30 + CH4 10 of 46
    assert rows.loc['B', 'n_starts'] == 2 and rows.loc['B', 'starts'] == 'start B'
    top = mr.largest_differences(tables.alternatives, (400.0, 240.0), 0.15)
    assert list(top.columns) == ['Species', 'A', 'B'] and len(top) <= 5
    assert top.Species.iloc[0] == 'Gas · carbon dioxide' and top.B.iloc[0] == pytest.approx(30 * 0.15)
    assert top.Species.iloc[1] == 'Solid · carbon'


def _texts(element):
    return [e.text for e in element.descendants() if isinstance(e, ui.label)]


def test_dialog_choices_and_comparison(tmp_path, config_dir):
    folder = stand_in_run(tmp_path / 'run')
    state = SavedState(folder, Settings.load(config_dir))
    holder = {}

    def page():
        holder['tab'] = ResultsTab(state)

    async def main():
        from nicegui.testing.user_simulation import user_simulation

        async with user_simulation(page) as user:
            await user.open('/')
            tab = holder['tab']
            with user:
                tab.select('saved:run')
            # banner: 2 of the 4 points
            await user.should_see('2 of 6 points ended at more than one result (different starting points gave '
                                  'different minima).')
            user.find('Compare them').click()
            dialog = tab.results_dialog
            await user.should_see(mr.TITLE)
            await user.should_see(mr.CAPTION)
            await user.should_see('400 °C · 240 bar')
            await user.should_see('500 °C · 240 bar')
            texts = _texts(dialog.detail)
            assert texts[:7] == list(mr.COLUMNS)
            assert 'Largest differences [mol / kg total feed]' in texts
            assert list(dialog.radios) == ['A', 'B'] and dialog.radios['A'].value == 'A'
            assert not dialog.checks['A'].enabled and dialog.checks['B'].value
            # choose B at 400 °C; at 500 °C show only B when comparing (untick C)
            with user:
                dialog.radios['B'].set_value('B')
            assert dialog.radios['A'].value is None and dialog.choices[(400.0, 240.0)] == 'B'
            user.find(marker='point-1').click()  # the list item '500 °C · 240 bar'
            assert list(dialog.radios) == ['A', 'B', 'C']
            with user:
                dialog.checks['C'].set_value(False)
            user.find('Apply choices').click()

            # choices.toml written and read back; tables use B at 400 °C
            assert read_choices(folder) == {(400.0, 240.0): 'B'}
            assert 'T_C = 400.0' in (folder / 'choices.toml').read_text()
            row = tab.view_tables().summary.set_index('T_C').loc[400.0]
            assert row.result == 'B' and row.gibbs_J_per_mol_C == pytest.approx(-1e6 + 120.0)
            gas = tab.view_tables().species
            co2 = gas[(gas.T_C == 400) & (gas.phase == 'gas') & (gas.species == 'carbon dioxide')]
            assert float(co2.mol_per_kg_dry_feed.iloc[0]) == pytest.approx(30.0)
            assert any(r.get('result') == 'B' for r in tab.summary_table.rows)
            with user:  # reading the run from disk again keeps the choice
                tab.load()
            assert tab.view_tables().choices == {(400.0, 240.0): 'B'}

            # comparison mode
            user.find('Compare them').click()
            assert tab.results_dialog.radios['B'].value == 'B'  # the dialog shows the saved choice
            user.find(marker='point-1').click()
            assert not tab.results_dialog.checks['C'].value  # the ticks are kept
            user.find(kind=ui.button, content='Compare in plots').click()
            await user.should_see('Comparing results at 2 points: A (solid, lowest Gibbs energy) and B (dashed). '
                                  'Elsewhere there is only one result.')
            spec = tab.spec
            assert [(s.dash, s.label) for s in spec.line_styles] == [('solid', 'Result A (lowest Gibbs energy)'),
                                                                     ('dash', 'Result B')]
            gas = spec.panels[0]
            base = {s.name: s for s in gas.series if s.legend}
            assert base['carbon dioxide'].dash == 'solid'
            assert base['carbon dioxide'].y[1] == pytest.approx(20.0)  # A at 400 °C, not the choice
            extra = {s.name: s for s in gas.series if s.lines and not s.legend}
            co2_b = extra['carbon dioxide (B)']
            assert co2_b.dash == 'dash'
            # x = 300, 400, 450 (failed), 500, 600, 700; B exists at 400 and 500 and is drawn from the
            # previous to the next point around them: 300 to 600, not at 700
            assert gas.series[0].x == [300.0, 400.0, 450.0, 500.0, 600.0, 700.0]
            assert co2_b.y[0] == pytest.approx(10.0) and co2_b.y[1] == pytest.approx(30.0)
            assert math.isnan(co2_b.y[2]) and co2_b.y[3] == pytest.approx(25.0) and co2_b.y[4] == pytest.approx(25.0)
            assert math.isnan(co2_b.y[5])
            methane_b = extra['methane (B)']
            assert methane_b.y[3] == pytest.approx(30.0) and methane_b.y[4] == pytest.approx(15.0)
            assert math.isnan(methane_b.y[5])
            assert not any(name.endswith('(C)') for name in extra)  # C is unticked
            marks = [s for s in gas.series if s.marker == 'circle-open' and s.name == 'carbon dioxide (B)']
            assert marks[0].x == [400.0, 500.0]
            assert [(b.x0, b.x1) for b in gas.bands] == [(300.0, 600.0)]
            figure = tab.plot.figure
            names = [t.name for t in figure.data]
            assert 'Result A (lowest Gibbs energy)' in names and 'Result B' in names
            # exports carry the line-style legend: the dialog saves this spec
            user.find('Save plot…').click()
            assert tab.save_dialog.spec.line_styles == spec.line_styles
            with user:
                tab.save_dialog.dialog.close()
            # Change reopens the dialog; Stop comparing goes back to one line per species
            user.find('Change').click()
            assert tab.results_dialog.dialog.value
            with user:
                tab.results_dialog.dialog.close()
            user.find('Stop comparing').click()
            await user.should_see('2 of 6 points ended at more than one result')
            assert tab.spec.line_styles == [] and not any(s.name.endswith('(B)') for s in tab.spec.panels[0].series)

            # back to A everywhere: choices.toml is removed
            user.find('Compare them').click()
            with user:
                tab.results_dialog.radios['A'].set_value('A')
            user.find('Apply choices').click()
            assert not (folder / 'choices.toml').exists() and tab.view_tables().choices == {}

    asyncio.run(main())
    assert AMOUNTS[400][('gas', 'carbon dioxide')] == 20.0  # the shared fixture data are unchanged
