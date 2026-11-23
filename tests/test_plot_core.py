"""Results-tab quantities, plot specifications and exports (WP5), without the GUI or Julia."""
import math
import re
import struct
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from htg import plotbuild as pb
from htg import plotrender as pr
from htg import plotspec as ps
from htg import quantities as q

FEED_C = 40.0  # mol feed carbon per kg dry feed
MW = {'water': 18.0, 'carbon dioxide': 44.0, 'hydrogen': 2.0, 'methane': 16.0, 'carbon': 12.0}
N_C = {'water': 0, 'carbon dioxide': 1, 'hydrogen': 0, 'methane': 1, 'carbon': 1}
# (T, phase, species): mol / kg dry feed; 450 °C failed; no organic liquid anywhere
AMOUNTS = {
    300: {('gas', 'water'): 50.0, ('gas', 'carbon dioxide'): 10.0, ('gas', 'hydrogen'): 2.0,
          ('aqueous', 'water'): 400.0, ('aqueous', 'carbon dioxide'): 5.0, ('solid', 'carbon'): 20.0},
    400: {('gas', 'water'): 300.0, ('gas', 'carbon dioxide'): 20.0, ('gas', 'hydrogen'): 5.0,
          ('gas', 'methane'): 10.0, ('aqueous', 'water'): 100.0, ('aqueous', 'carbon dioxide'): 2.0,
          ('solid', 'carbon'): 8.0},
    500: {('gas', 'water'): 450.0, ('gas', 'carbon dioxide'): 25.0, ('gas', 'hydrogen'): 15.0,
          ('gas', 'methane'): 15.0},
}


def make_tables(amounts=AMOUNTS, P=240.0, failed=(450,)):
    rows, phase_rows, summary = [], [], []
    for T, items in amounts.items():
        carbon = {p: 0.0 for p in ('gas', 'aqueous', 'organic_liquid', 'solid')}
        for (phase, name), mol in items.items():
            rows.append({'T_C': float(T), 'P_bar': P, 'phase': phase,
                         'group': 'aqueous neutrals' if phase == 'aqueous' else phase, 'species': name,
                         'formula': '', 'mol_per_mol_C': mol / FEED_C, 'mol_per_kg_dry_feed': mol,
                         'g_per_kg_dry_feed': mol * MW[name], 'carbon_mol_per_kg_dry_feed': mol * N_C[name],
                         'molality_mol_per_kg_water': np.nan})
            carbon[phase] += mol * N_C[name]
        total_C = sum(carbon.values())
        for phase, c in carbon.items():
            phase_rows.append({'T_C': float(T), 'P_bar': P, 'phase': phase, 'mol_per_mol_C': 0.0,
                               'g_per_kg_dry_feed': 0.0, 'carbon_percent': 100 * c / total_C,
                               'gibbs_J_per_mol_C': 0.0})
        summary.append({'T_C': float(T), 'P_bar': P, 'status': 'solved', 'gibbs_J_per_mol_C': -1e6,
                        'gas_carbon_percent': 100 * carbon['gas'] / total_C})
    for T in failed:
        summary.append({'T_C': float(T), 'P_bar': P, 'status': 'failed', 'gibbs_J_per_mol_C': np.nan,
                        'gas_carbon_percent': np.nan})
    summary = pd.DataFrame(summary).sort_values('T_C').reset_index(drop=True)
    return SimpleNamespace(summary=summary, species=pd.DataFrame(rows), phases=pd.DataFrame(phase_rows))


@pytest.fixture
def tables():
    return make_tables()


def value(long, T, name, phase=None):
    part = long[(long.T_C == T) & (long.species == name)]
    if phase is not None:
        part = part[part.phase == phase]
    return float(part.value.sum())


# --- quantities -------------------------------------------------------------------------------

def test_unit_bases(tables):
    s, a = 0.2, 0.1
    dry = q.amounts(tables.species, 'mol', 'dry')
    total = q.amounts(tables.species, 'mol', 'total', s, a)
    daf = q.amounts(tables.species, 'mol', 'daf', s, a)
    assert value(dry, 400, 'carbon dioxide', 'gas') == 20.0
    assert value(total, 400, 'carbon dioxide', 'gas') == pytest.approx(20.0 * s)
    assert value(daf, 400, 'carbon dioxide', 'gas') == pytest.approx(20.0 / (1 - a))
    grams = q.amounts(tables.species, 'g', 'total', s, a)
    assert value(grams, 400, 'carbon dioxide', 'gas') == pytest.approx(20.0 * 44.0 * s)
    with pytest.raises(ValueError, match='no total mass'):
        q.amounts(tables.species, 'mol', 'total', None, a)
    assert q.available_bases(None)['total'] == q.NO_TOTAL_MASS
    assert q.available_bases(0.2) == {'total': None, 'dry': None, 'daf': None}
    assert q.amount_unit('g', 'daf') == 'g / kg dry ash-free feed'
    assert q.amount_unit('mol', 'total') == 'mol / kg total feed'


def test_all_phases_together(tables):
    total = q.total_amounts(tables.species, 'mol', 'total', 0.5)
    assert value(total, 300, 'carbon dioxide') == pytest.approx((10.0 + 5.0) * 0.5)
    assert value(total, 300, 'water') == pytest.approx(450.0 * 0.5)
    assert len(total[(total.T_C == 300) & (total.species == 'carbon dioxide')]) == 1


@pytest.mark.parametrize('by', ['mol', 'mass'])
@pytest.mark.parametrize('dry', [False, True])
def test_gas_composition(tables, by, dry):
    comp = q.gas_composition(tables.species, by, exclude_water=dry)
    at400 = comp[comp.T_C == 400].set_index('species').value
    assert at400.sum() == pytest.approx(100.0)
    mol = {'water': 300.0, 'carbon dioxide': 20.0, 'hydrogen': 5.0, 'methane': 10.0}
    if dry:
        mol.pop('water')
        assert 'water' not in at400
    weight = {k: v * (MW[k] if by == 'mass' else 1.0) for k, v in mol.items()}
    for name, w in weight.items():
        assert at400[name] == pytest.approx(100 * w / sum(weight.values()))


def test_cge_is_gas_carbon_percent(tables):
    cge = q.carbon_gasification_efficiency(tables.summary)
    solved = tables.summary[tables.summary.status == 'solved']
    assert cge.value.tolist() == pytest.approx(solved.gas_carbon_percent.tolist())
    # only gas-phase carbon: at 300 °C gas 10, aqueous 5, solid 20 mol C
    assert float(cge[cge.T_C == 300].value.iloc[0]) == pytest.approx(100 * 10 / 35)
    spec = pb.build_plot('cge', tables)
    assert spec.panels[0].series[0].y[:2] == pytest.approx([100 * 10 / 35, 100 * 30 / 40])


def test_profile_gaps_and_zeros(tables):
    long = q.amounts(tables.species)
    table = q.profile(long[long.phase == 'gas'], tables.summary, 'T_C', 240.0, ['methane'])
    assert table.index.tolist() == [300.0, 400.0, 450.0, 500.0]
    assert table['methane'].iloc[0] == 0.0  # solved, no methane
    assert math.isnan(table['methane'].iloc[2])  # failed point


def test_measured_amounts(tables):
    measured = pd.DataFrame({'T_C': [400.0], 'gas:carbon dioxide': [0.5]})
    long = q.measured_amounts(measured, 'mol_per_mol_C', tables.species, 'g', 'total', 0.5)
    assert float(long.value.iloc[0]) == pytest.approx(0.5 * FEED_C * 44.0 * 0.5)


# --- plot specification -----------------------------------------------------------------------

def test_phases_spec_flags_empty_organic_liquid(tables):
    opts = pb.PlotOptions(basis='total', s=0.15)
    spec = pb.build_plot('phases', tables, opts)
    assert spec.title == 'Amount in each phase · mol / kg total feed · 240 bar'
    assert [p.title for p in spec.panels] == ['Gas', 'Aqueous', 'Solid', 'Organic liquid']
    assert (spec.rows, spec.cols) == (2, 2)
    assert [p.empty for p in spec.panels] == [False, False, False, True]
    assert spec.panels[3].empty_message == 'No organic liquid at any point'
    gas = {s.name: s for s in spec.panels[0].series}
    assert 'water' not in gas  # water off by default
    assert gas['carbon dioxide'].y == pytest.approx([1.5, 3.0, math.nan, 3.75], nan_ok=True)
    keys = [o.key for o in ps.layout_options(spec)]
    assert keys == ['one', 'grid', 'row', 'present']
    present = ps.layout_options(spec)[3]
    assert present.caption == 'Gas, Aqueous, Solid (no organic liquid at any point)'
    assert present.panels == [0, 1, 2]
    assert ps.summary_line(spec, 'present', 190, 65, 'pdf') == '3 panels · 190 × 65 mm · PDF'
    assert ps.summary_line(spec, 'one', 85, 65, 'png') == '1 panel · 85 × 65 mm · PNG'


def test_single_panel_plots_offer_one_layout(tables):
    for kind in ('gas', 'cge', 'carbon', 'total'):
        spec = pb.build_plot(kind, tables)
        assert [o.key for o in ps.layout_options(spec)] == ['one']
    full = make_tables({T: {**items, ('organic_liquid', 'methane'): 1.0} for T, items in AMOUNTS.items()})
    spec = pb.build_plot('phases', full)
    assert [o.key for o in ps.layout_options(spec)] == ['one', 'grid', 'row']


def test_colour_repetition_gives_dashes():
    n = len(ps.PALETTES['Default'])
    first = ps.resolve_style(ps.Series([], [], 'a', colour=0), 'Default')
    again = ps.resolve_style(ps.Series([], [], 'b', colour=n), 'Default')
    third = ps.resolve_style(ps.Series([], [], 'c', colour=2 * n), 'Default')
    assert first == ('#1f77b4', 'solid')
    assert again == ('#1f77b4', 'dash')
    assert third == ('#1f77b4', 'dot')
    assert ps.resolve_style(ps.Series([], [], 'd', colour=n), 'Set3') == (ps.PALETTES['Set3'][n], 'solid')
    assert ps.resolve_style(ps.Series([], [], 'e', colour=n, dash='dashdot'), 'Default')[1] == 'dashdot'
    # a run with more species than colours
    many = {300: {('gas', f's{i}'): float(20 - i) for i in range(10)}}
    MW.update({f's{i}': 10.0 for i in range(10)})
    N_C.update({f's{i}': 1 for i in range(10)})
    spec = pb.build_plot('total', make_tables(many, failed=()), pb.PlotOptions(species=[f's{i}' for i in range(10)]))
    dashes = [ps.resolve_style(s, 'Default')[1] for s in spec.panels[0].series]
    assert dashes == ['solid'] * 8 + ['dash'] * 2


def test_colours_stable_across_plots(tables):
    opts = pb.PlotOptions(species=['carbon dioxide', 'hydrogen', 'methane', 'carbon'])
    colours = {}
    for kind in ('total', 'phases', 'gas'):
        for panel in pb.build_plot(kind, tables, opts).panels:
            for s in panel.series:
                colours.setdefault(s.key, set()).add(s.colour)
    assert all(len(c) == 1 for c in colours.values())
    assert pb.run_colours(tables.species)['water'] == max(pb.run_colours(tables.species).values())
    assert list(ps.PALETTES) == ['Default', 'Set1', 'Set2', 'Set3', 'Dark2', 'Accent', 'Paired', 'Pastel1',
                                 'Pastel2']
    assert ps.LINE_WIDTHS == (1, 1.5, 2, 3, 4)


def test_largest_and_water(tables):
    assert pb.largest_species('total', tables, pb.PlotOptions(), n=2) == ['carbon dioxide', 'carbon']
    spec = pb.build_plot('total', tables, pb.PlotOptions(include_water=True))
    assert 'water' in [s.name for s in spec.panels[0].series]
    spec = pb.build_plot('gas', tables, pb.PlotOptions(include_water=True, exclude_water=True))
    assert 'water' not in [s.name for s in spec.panels[0].series]


def test_units_and_cutoff(tables):
    opts = pb.PlotOptions(cutoff_C=390.0, t_unit='K', t_convert=lambda T: T + 273.15, p_unit='MPa',
                          p_convert=lambda P: P / 10)
    spec = pb.build_plot('phases', tables, opts)
    assert spec.title.endswith(' · 24 MPa')
    gas, aqueous = spec.panels[0], spec.panels[1]
    assert gas.x_label == 'Temperature [K]'
    assert gas.series[0].x[0] == pytest.approx(573.15)
    assert gas.vlines == []
    assert aqueous.vlines[0].x == pytest.approx(663.15)
    assert aqueous.vlines[0].label == 'aqueous removed ≥ 663.15 K'
    assert aqueous.vlines[0].dash == 'dash'
    # pressure axis: no cut-off line, x are pressures at the fixed temperature
    grid = make_tables()
    other = make_tables(P=300.0)
    both = SimpleNamespace(**{k: pd.concat([getattr(grid, k), getattr(other, k)], ignore_index=True)
                              for k in ('summary', 'species', 'phases')})
    opts.x_axis, opts.fixed = 'P', 400.0
    spec = pb.build_plot('total', both, opts)
    assert spec.panels[0].x_label == 'Pressure [MPa]'
    assert spec.panels[0].series[0].x == pytest.approx([24.0, 30.0])
    assert spec.panels[0].vlines == []
    assert spec.title.endswith(' · 673.15 K')
    assert pb.map_available(both.summary)
    heat = pb.build_plot('map', both, pb.PlotOptions(map_item='species:gas:methane')).panels[0].heatmap
    assert heat.x == ['300', '400', '450', '500'] and heat.y == ['240', '300']
    assert heat.z[0][0] == 0.0 and math.isnan(heat.z[0][2]) and heat.z[0][1] == 10.0


def test_log_scale_drops_zeros(tables):
    spec = pb.build_plot('total', tables, pb.PlotOptions(log=True, species=['methane']))
    assert spec.panels[0].log
    assert math.isnan(spec.panels[0].series[0].y[0])


def test_measured_series(tables):
    measured = pb.Measured(pd.DataFrame({'T_C': [300.0, 400.0], 'gas:carbon dioxide': [9.0, 21.0]}),
                           'mol_per_kg_dry_feed', 'Experiment')
    spec = pb.build_plot('phases', tables, pb.PlotOptions(measured=measured))
    marks = [s for s in spec.panels[0].series if s.measured]
    assert marks[0].y == [9.0, 21.0] and marks[0].marker == 'square-open' and not marks[0].lines
    assert marks[-1].name == 'Experiment' and marks[-1].legend
    assert not any(s.measured for s in spec.panels[1].series)


def test_comparison_spec(tables):
    opts = pb.PlotOptions(species=['carbon dioxide', 'methane'])
    base = pb.build_plot('phases', tables, opts)
    changed = {T: dict(items) for T, items in AMOUNTS.items()}
    changed[400][('gas', 'carbon dioxide')] = 30.0
    other = pb.build_plot('phases', make_tables(changed), opts)
    flagged = pb.flagged_x([(400.0, 240.0)], tables, opts)
    spec = ps.compare_results(base, {'B': other}, flagged, descriptions={'A': 'lowest Gibbs energy'})
    gas = spec.panels[0]
    extra = [s for s in gas.series if s.name == 'carbon dioxide (B)' and s.lines]
    assert len(extra) == 1 and extra[0].dash == 'dash' and not extra[0].legend
    # drawn from the previous to the next point only: 300, 400, 450 (failed, gap); 500 left out
    assert extra[0].y[:2] == [10.0, 30.0] and math.isnan(extra[0].y[3])
    marks = [s for s in gas.series if s.marker == 'circle-open']
    assert marks[0].x == [400.0] and marks[0].y == [30.0]
    assert all(s.dash == 'solid' for s in gas.series if s.legend and s.lines)
    assert [(b.x0, b.x1) for b in gas.bands] == [(300.0, 450.0)]
    assert [(s.dash, s.label) for s in spec.line_styles] == [('solid', 'Result A (lowest Gibbs energy)'),
                                                             ('dash', 'Result B')]
    assert pr.line_style_text(spec.line_styles[1]) == '- - Result B'
    assert spec.panels[3].empty and base.panels[0].bands == []


def test_no_points():
    empty = SimpleNamespace(summary=pd.DataFrame(columns=['T_C', 'P_bar', 'status']),
                            species=pd.DataFrame(), phases=pd.DataFrame())
    for kind in pb.PLOTS:
        spec = pb.build_plot(kind, empty)
        assert spec.message == pb.NO_POINTS
        pr.render_plotly(spec)


# --- rendering --------------------------------------------------------------------------------

def test_render_plotly(tables):
    spec = pb.build_plot('phases', tables, pb.PlotOptions(cutoff_C=390.0, log=True))
    fig = pr.render_plotly(spec, 'Dark2', 3)
    traces = [t for t in fig.data if t.name == 'carbon dioxide']
    assert len(traces) == 2 and sum(bool(t.showlegend) for t in traces) == 1
    assert traces[0].line.color == ps.PALETTES['Dark2'][pb.run_colours(tables.species)['carbon dioxide']]
    assert traces[0].line.width == 3
    texts = [a.text for a in fig.layout.annotations]
    assert 'No organic liquid at any point' in texts and 'aqueous removed ≥ 390 °C' in texts
    assert fig.layout.yaxis.type == 'log'
    single = pr.render_plotly(pb.build_plot('cge', tables))
    assert single.layout.title.text.startswith('Carbon gasification efficiency · % · 240 bar')
    pr.render_plotly(pb.build_plot('map', tables))


def _png_size(path):
    data = path.read_bytes()[:24]
    return struct.unpack('>II', data[16:24])


def _svg_size_mm(path):
    head = path.read_text()[:2000]
    w = float(re.search(r'width="([\d.]+)pt"', head).group(1))
    h = float(re.search(r'height="([\d.]+)pt"', head).group(1))
    return w / 72 * 25.4, h / 72 * 25.4


def _pdf_size_mm(path):
    box = re.search(rb'/MediaBox\s*\[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]', path.read_bytes())
    x0, y0, x1, y1 = map(float, box.groups())
    return (x1 - x0) / 72 * 25.4, (y1 - y0) / 72 * 25.4


@pytest.mark.parametrize('layout', ['one', 'grid', 'row', 'present'])
@pytest.mark.parametrize('fmt', ['PNG', 'SVG', 'PDF'])
def test_export_size(tables, tmp_path, layout, fmt):
    spec = pb.build_plot('phases', tables, pb.PlotOptions(cutoff_C=390.0))
    width, height = ps.DEFAULT_SIZES_MM[layout]
    if layout == 'one':
        width, height = 91.5, 70.3
    dpi = 120
    path = tmp_path / f'plot.{fmt.lower()}'
    pr.render_matplotlib(spec, path, layout, width, height, fmt, dpi=dpi, panel=1, theme='Paired', line_width=1.5)
    if fmt == 'PNG':
        w, h = _png_size(path)
        assert abs(w - width / 25.4 * dpi) <= 1 and abs(h - height / 25.4 * dpi) <= 1
    else:
        w, h = (_svg_size_mm if fmt == 'SVG' else _pdf_size_mm)(path)
        assert abs(w - width) <= 0.5 and abs(h - height) <= 0.5
        assert b'<image' not in path.read_bytes()[:200000] if fmt == 'SVG' else True  # vector


@pytest.mark.parametrize('kind', ['total', 'gas', 'cge', 'carbon', 'map'])
def test_export_single_panel_plots(tables, tmp_path, kind):
    spec = pb.build_plot(kind, tables)
    path = tmp_path / 'plot.pdf'
    pr.render_matplotlib(spec, path, 'one', shared_legend=False)
    assert _pdf_size_mm(path) == pytest.approx(ps.DEFAULT_SIZES_MM['one'], abs=0.5)


def test_export_png_300_dpi_and_file_object(tables, tmp_path):
    import io
    spec = ps.compare_results(pb.build_plot('total', tables), {'B': pb.build_plot('total', tables)}, [400.0])
    path = tmp_path / 'plot.png'
    pr.render_matplotlib(spec, path, 'one', 85, 65)
    assert _png_size(path) == pytest.approx((85 / 25.4 * 300, 65 / 25.4 * 300), abs=1)
    buffer = io.BytesIO()
    pr.render_matplotlib(spec, buffer, 'one', format='svg')
    assert buffer.getvalue().lstrip().startswith(b'<?xml')
