"""Plot specification: plots as plain data, independent of any plotting library.

A PlotSpec holds everything needed to draw a plot (panels, series, labels, marker
lines, shaded ranges). It is built from result tables by htg.plotbuild and drawn
by htg.plotrender, on screen with Plotly and for export with Matplotlib, from the
same object. Nothing here imports NiceGUI, Plotly or Matplotlib.

Colours and dashes are resolved only when drawing (resolve_style), so changing the
colour theme or the line width does not rebuild the spec:

- Series.colour is a palette index (stable per species name, see assign_colours) or
  an explicit '#rrggbb'.
- Series.dash None means automatic: solid, or dashed then dotted when the palette
  has fewer colours than the index needs. An explicit dash ('solid', 'dash', 'dot',
  'dashdot') is kept as given; comparison mode uses this for results A, B, C, D.
"""
import copy
import math
from dataclasses import dataclass, field

# --- look ---------------------------------------------------------------------------------------

PALETTES = {
    'Default': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#4d4d4d', '#e377c2'],
    'Set1': ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf',
             '#999999'],
    'Set2': ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3'],
    'Set3': ['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5',
             '#d9d9d9', '#bc80bd', '#ccebc5', '#ffed6f'],
    'Dark2': ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666'],
    'Accent': ['#7fc97f', '#beaed4', '#fdc086', '#ffff99', '#386cb0', '#f0027f', '#bf5b17', '#666666'],
    'Paired': ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00',
               '#cab2d6', '#6a3d9a', '#ffff99', '#b15928'],
    'Pastel1': ['#fbb4ae', '#b3cde3', '#ccebc5', '#decbe4', '#fed9a6', '#ffffcc', '#e5d8bd', '#fddaec',
                '#f2f2f2'],
    'Pastel2': ['#b3e2cd', '#fdcdac', '#cbd5e8', '#f4cae4', '#e6f5c9', '#fff2ae', '#f1e2cc', '#cccccc'],
}
THEME_LABELS = {name: f'{name} (light)' if name in ('Set3', 'Pastel1', 'Pastel2') else name
                for name in PALETTES}
DEFAULT_THEME = 'Default'
LINE_WIDTHS = (1, 1.5, 2, 3, 4)  # px on screen; exports use 0.75 pt per px
DEFAULT_LINE_WIDTH = 2
DASHES = ('solid', 'dash', 'dot', 'dashdot')
REPEAT_DASHES = ('solid', 'dash', 'dot')  # first, second and third round through a palette
RESULT_DASHES = ('solid', 'dash', 'dot', 'dashdot')  # results A, B, C, D in comparison mode
DASH_SYMBOLS = {'solid': '—', 'dash': '- -', 'dot': '···', 'dashdot': '-·-'}
NEUTRAL = '#333333'  # line-style legend and single neutral lines


def palette(theme):
    return PALETTES.get(theme, PALETTES[DEFAULT_THEME])


def assign_colours(names, start=0):
    """{name: palette index} in the given order (duplicates ignored).

    Use one assignment per run for every plot and panel so a species keeps its colour.
    """
    out = {}
    for name in names:
        if name not in out:
            out[name] = start + len(out)
    return out


def resolve_style(series, theme=DEFAULT_THEME):
    """(colour '#rrggbb', dash) of a series for the theme.

    Palette indices beyond the palette length wrap around; the second round is dashed
    and the third dotted (unless the series has an explicit dash).
    """
    colours = palette(theme)
    if isinstance(series.colour, str):
        return series.colour, series.dash or 'solid'
    index = int(series.colour)
    colour = colours[index % len(colours)]
    dash = series.dash or REPEAT_DASHES[(index // len(colours)) % len(REPEAT_DASHES)]
    return colour, dash


# --- the specification --------------------------------------------------------------------------

@dataclass
class Series:
    """One line (or set of markers) in a panel.

    x, y      lists of floats; NaN in y is a gap (failed point)
    name      legend text
    colour    palette index (int) or explicit '#rrggbb'
    key       legend group; series with the same key share one legend entry (species name,
              phase name, ...)
    dash      None (automatic, see resolve_style) or one of DASHES
    marker    None, 'circle', 'circle-open' or 'square-open'
    lines     False draws markers only
    measured  measured (reference) data rather than a calculation
    legend    False keeps the series out of the legend
    """
    x: list
    y: list
    name: str
    colour: int | str = 0
    key: str = ''
    dash: str | None = None
    marker: str | None = 'circle'
    lines: bool = True
    measured: bool = False
    legend: bool = True


@dataclass
class VLine:
    """Vertical marker line, e.g. the aqueous cut-off (x in the panel's x unit)."""
    x: float
    label: str = ''
    dash: str = 'dash'


@dataclass
class Band:
    """Lightly shaded x range, e.g. points with more than one result."""
    x0: float
    x1: float
    label: str = ''


@dataclass
class Heatmap:
    """Colour map over a grid; z[i][j] belongs to y[i], x[j]; NaN = no value."""
    x: list
    y: list
    z: list
    colourbar: str = ''
    colourscale: str = 'Viridis'


@dataclass
class LineStyle:
    """Entry of the line-style legend, e.g. LineStyle('dash', 'Result B')."""
    dash: str
    label: str


@dataclass
class Panel:
    """One set of axes.

    empty_message  when set, the panel shows this text instead of a plot
    key            identifies the panel (phase name for "Amount in each phase")
    """
    title: str
    x_label: str = ''
    y_label: str = ''
    series: list = field(default_factory=list)
    log: bool = False
    empty_message: str | None = None
    vlines: list = field(default_factory=list)
    bands: list = field(default_factory=list)
    heatmap: Heatmap | None = None
    key: str = ''

    @property
    def empty(self):
        return self.empty_message is not None


@dataclass
class PlotSpec:
    """A whole plot: title, panels on a rows x cols grid, optional line-style legend.

    kind     the plot type key (see htg.plotbuild.PLOTS)
    message  shown instead of the plot when there is nothing to draw (e.g. no points yet)
    """
    title: str
    panels: list
    kind: str = ''
    rows: int = 1
    cols: int = 1
    line_styles: list = field(default_factory=list)
    message: str = ''


# --- export layouts -----------------------------------------------------------------------------

@dataclass
class LayoutOption:
    key: str
    label: str
    caption: str
    panels: list  # indices into spec.panels ("one": all choices for the Panel select)
    size_mm: tuple


LAYOUTS = {'one': 'One panel', 'grid': 'All four, 2 × 2', 'row': 'All four, in a row',
           'present': 'Phases present, in a row'}
LAYOUT_CAPTIONS = {'one': '', 'grid': 'As shown on the Results tab', 'row': 'One row of four panels'}
DEFAULT_SIZES_MM = {'one': (85, 65), 'grid': (170, 130), 'row': (250, 65), 'present': (190, 65)}
FORMATS = ('PNG', 'SVG', 'PDF')
DEFAULT_DPI = 300


def layout_options(spec):
    """Export layouts offered for a spec, in display order.

    Single-panel plots only offer "One panel". "Phases present, in a row" is offered
    only when at least one panel is empty and at least one is not.
    """
    indices = list(range(len(spec.panels)))
    options = [LayoutOption('one', LAYOUTS['one'], '', indices, DEFAULT_SIZES_MM['one'])]
    if len(spec.panels) < 4:
        return options
    options.append(LayoutOption('grid', LAYOUTS['grid'], LAYOUT_CAPTIONS['grid'], indices,
                                DEFAULT_SIZES_MM['grid']))
    options.append(LayoutOption('row', LAYOUTS['row'], LAYOUT_CAPTIONS['row'], indices,
                                DEFAULT_SIZES_MM['row']))
    present = [i for i, p in enumerate(spec.panels) if not p.empty]
    if 0 < len(present) < len(spec.panels):
        options.append(LayoutOption('present', LAYOUTS['present'], present_caption(spec), present,
                                    DEFAULT_SIZES_MM['present']))
    return options


def present_caption(spec):
    """'Gas, Aqueous, Solid (no organic liquid at any point)'."""
    present = [p.title for p in spec.panels if not p.empty]
    missing = [p.title.lower() for p in spec.panels if p.empty]
    if not missing:
        return ', '.join(present)
    gone = missing[0] if len(missing) == 1 else ', '.join(missing[:-1]) + ' or ' + missing[-1]
    return f"{', '.join(present)} (no {gone} at any point)"


def layout_panels(spec, layout='one', panel=0):
    """(panel indices, rows, cols) of an export layout."""
    if layout == 'one' or len(spec.panels) == 1:
        return [min(max(int(panel), 0), len(spec.panels) - 1)], 1, 1
    if layout == 'grid':
        cols = 2
        return list(range(len(spec.panels))), math.ceil(len(spec.panels) / cols), cols
    if layout == 'row':
        return list(range(len(spec.panels))), 1, len(spec.panels)
    if layout == 'present':
        present = [i for i, p in enumerate(spec.panels) if not p.empty] or [0]
        return present, 1, len(present)
    raise ValueError(f'unknown layout {layout!r}')


def summary_line(spec, layout, width_mm, height_mm, fmt, panel=0):
    """'3 panels · 190 × 65 mm · PDF'."""
    n = len(layout_panels(spec, layout, panel)[0])
    return f"{n} panel{'s' if n != 1 else ''} · {width_mm:g} × {height_mm:g} mm · {fmt.upper()}"


# --- comparison of results (WP6) ----------------------------------------------------------------

def result_label(label, description=''):
    text = f'Result {label}'
    return f'{text} ({description})' if description else text


def compare_results(base, others, flagged, base_label='A', descriptions=None):
    """Comparison-mode spec: result A everywhere, other results only around flagged points.

    base         PlotSpec built from the tables of result A (the chosen result)
    others       {label: PlotSpec} built the same way (same options) from the tables with the
                 flagged points replaced by that result (htg.quantities.alternative_tables)
    flagged      x values (in the plot's x unit) of the points with more than one result
    descriptions {label: text} for the line-style legend, e.g. {'A': 'lowest Gibbs energy'}

    Result A is drawn solid; B dashed, C dotted, D dash-dot. Each other result is drawn
    only from the previous to the next point around each flagged point, with open circles
    at the flagged points. The flagged ranges are shaded and a line-style legend is added.
    Series are matched by panel and (key, name).
    """
    descriptions = descriptions or {}
    spec = copy.deepcopy(base)
    flagged = [float(x) for x in flagged]
    labels = [base_label, *others]
    spec.line_styles = [LineStyle(RESULT_DASHES[i % len(RESULT_DASHES)],
                                  result_label(label, descriptions.get(label, '')))
                        for i, label in enumerate(labels)]
    for panel in spec.panels:
        for series in panel.series:
            if series.lines and not series.measured:
                series.dash = RESULT_DASHES[0]
    for i, (label, other) in enumerate(others.items(), start=1):
        dash = RESULT_DASHES[i % len(RESULT_DASHES)]
        for panel, alt_panel in zip(spec.panels, other.panels):
            if panel.heatmap is not None or panel.empty:
                continue
            for series in alt_panel.series:
                if not series.lines or series.measured:
                    continue
                y_window, x_marks, y_marks = _window(series.x, series.y, flagged)
                if not any(_finite(v) for v in y_window):
                    continue
                panel.series.append(Series(list(series.x), y_window, f'{series.name} ({label})',
                                           series.colour, series.key, dash, None, True, False, False))
                panel.series.append(Series(x_marks, y_marks, f'{series.name} ({label})', series.colour,
                                           series.key, dash, 'circle-open', False, False, False))
    for panel in spec.panels:
        if panel.heatmap is None and not panel.empty and panel.series:
            panel.bands.extend(_bands(panel.series[0].x, flagged))
    return spec


def _finite(v):
    try:
        return math.isfinite(float(v))
    except (TypeError, ValueError):
        return False


def _near(a, b):
    return abs(a - b) <= 1e-9 * max(1.0, abs(a), abs(b))


def _window(x, y, flagged):
    """y kept from the previous to the next point around each flagged x, NaN elsewhere."""
    keep = [False] * len(x)
    marks_x, marks_y = [], []
    for j, xj in enumerate(x):
        if any(_near(float(xj), f) for f in flagged):
            for k in (j - 1, j, j + 1):
                if 0 <= k < len(x):
                    keep[k] = True
            marks_x.append(xj)
            marks_y.append(y[j])
    return [v if k else math.nan for v, k in zip(y, keep)], marks_x, marks_y


def _bands(x, flagged):
    """Shaded ranges from the previous to the next point around flagged points, merged."""
    x = [float(v) for v in x]
    order = sorted(range(len(x)), key=lambda j: x[j])
    xs = [x[j] for j in order]
    ranges = []
    for j, xj in enumerate(xs):
        if any(_near(xj, f) for f in flagged):
            lo, hi = xs[max(j - 1, 0)], xs[min(j + 1, len(xs) - 1)]
            if ranges and lo <= ranges[-1][1]:
                ranges[-1][1] = max(ranges[-1][1], hi)
            else:
                ranges.append([lo, hi])
    return [Band(lo, hi) for lo, hi in ranges]
