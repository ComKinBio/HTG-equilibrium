"""Draw a PlotSpec: on screen with Plotly, for export with Matplotlib.

Both renderers take the colour theme (a key of htg.plotspec.PALETTES) and the
line width in px (one of LINE_WIDTHS). Exports use 0.75 pt per px, fonts of
8 pt (9 pt titles), the exact page size requested and vector output for SVG
and PDF.
"""
import math

from .plotspec import (
    DASH_SYMBOLS,
    DEFAULT_DPI,
    DEFAULT_LINE_WIDTH,
    DEFAULT_SIZES_MM,
    DEFAULT_THEME,
    NEUTRAL,
    layout_panels,
    resolve_style,
)

BAND_COLOUR = '#9e9e9e'
PT_PER_PX = 0.75


def _finite(v):
    try:
        return math.isfinite(float(v))
    except (TypeError, ValueError):
        return False


# --- Plotly -------------------------------------------------------------------------------------

_PLOTLY_MARKERS = {'circle': 'circle', 'circle-open': 'circle-open', 'square-open': 'square-open'}


def render_plotly(spec, theme=DEFAULT_THEME, line_width=DEFAULT_LINE_WIDTH):
    """Plotly Figure of a spec: one subplot per panel; empty panels show their message."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    rows, cols = max(spec.rows, 1), max(spec.cols, 1)
    multi = len(spec.panels) > 1
    if multi:
        fig = make_subplots(rows=rows, cols=cols, subplot_titles=[p.title for p in spec.panels],
                            horizontal_spacing=0.09, vertical_spacing=0.14)
    else:
        fig = go.Figure()
    shown = set()
    for i, panel in enumerate(spec.panels):
        where = {'row': i // cols + 1, 'col': i % cols + 1} if multi else {}
        if multi:
            sub = fig.get_subplot(**where)
            xref, yref = sub.xaxis.plotly_name.replace('axis', ''), sub.yaxis.plotly_name.replace('axis', '')
        else:
            xref, yref = 'x', 'y'
        if panel.empty:
            fig.update_xaxes(visible=False, **where)
            fig.update_yaxes(visible=False, **where)
            fig.add_annotation(text=panel.empty_message, x=0.5, y=0.5, xref=f'{xref} domain',
                               yref=f'{yref} domain', showarrow=False, font={'size': 13, 'color': '#757575'})
            continue
        if panel.heatmap is not None:
            h = panel.heatmap
            fig.add_trace(go.Heatmap(z=h.z, x=h.x, y=h.y, colorscale=h.colourscale, colorbar={'title': h.colourbar},
                                     hovertemplate='%{x}, %{y}<br>%{z:.4g}<extra></extra>'), **where)
            fig.update_xaxes(type='category', **where)
            fig.update_yaxes(type='category', **where)
        for band in panel.bands:
            fig.add_vrect(x0=band.x0, x1=band.x1, fillcolor=BAND_COLOUR, opacity=0.15, line_width=0,
                          layer='below', **where)
        for s in panel.series:
            colour, dash = resolve_style(s, theme)
            mode = 'lines+markers' if s.lines and s.marker else 'lines' if s.lines else 'markers'
            key = s.key or s.name
            show = s.legend and key not in shown
            if show:
                shown.add(key)
            marker = {'color': colour, 'size': 8 if s.measured or s.marker == 'circle-open' else 5}
            if s.marker:
                marker['symbol'] = _PLOTLY_MARKERS.get(s.marker, s.marker)
            fig.add_trace(go.Scatter(x=list(s.x), y=list(s.y), name=s.name, mode=mode, legendgroup=key,
                                     showlegend=show, connectgaps=False,
                                     line={'color': colour, 'dash': dash, 'width': line_width}, marker=marker,
                                     hovertemplate=f'{s.name}<br>%{{x:.6g}}<br>%{{y:.4g}}<extra></extra>'),
                          **where)
        for line in panel.vlines:
            fig.add_vline(x=line.x, line_dash=line.dash, line_color=NEUTRAL, line_width=1,
                          annotation_text=line.label, annotation_position='top left',
                          annotation_font={'size': 11, 'color': NEUTRAL}, **where)
        fig.update_xaxes(title_text=panel.x_label, **where)
        fig.update_yaxes(title_text=panel.y_label, **where)
        if panel.log:
            fig.update_yaxes(type='log', exponentformat='power', **where)
        if not multi and panel.title:
            fig.update_layout(title=f'{spec.title}<br><sup>{panel.title}</sup>')
    for style in spec.line_styles:
        fig.add_trace(go.Scatter(x=[None], y=[None], mode='lines', name=style.label, legendgroup='line styles',
                                 line={'color': NEUTRAL, 'dash': style.dash, 'width': line_width}))
    if spec.message:
        fig.add_annotation(text=spec.message, x=0.5, y=0.5, xref='paper', yref='paper', showarrow=False,
                           font={'size': 14, 'color': '#757575'})
    if multi or not spec.panels or not spec.panels[0].title:
        fig.update_layout(title=spec.title)
    fig.update_layout(template='plotly_white', legend={'font': {'size': 11}}, hovermode='closest',
                      margin={'l': 60, 'r': 20, 't': 70 if multi else 50, 'b': 50})
    return fig


# --- Matplotlib ---------------------------------------------------------------------------------

_MPL_DASHES = {'solid': '-', 'dash': '--', 'dot': ':', 'dashdot': '-.'}
_MPL_MARKERS = {'circle': ('o', True), 'circle-open': ('o', False), 'square-open': ('s', False)}
FONT_PT = 8
TITLE_PT = 9
MM_PER_INCH = 25.4


def render_matplotlib(spec, path, layout='one', width_mm=None, height_mm=None, format=None, dpi=DEFAULT_DPI,
                      shared_legend=True, theme=DEFAULT_THEME, line_width=DEFAULT_LINE_WIDTH, panel=0,
                      title=True):
    """Save a spec with Matplotlib at an exact page size.

    path       file path, or a binary file object (then ``format`` is required)
    layout     'one' (the panel given by ``panel``), 'grid' (2 × 2), 'row' (all in a row) or
               'present' (only the non-empty panels, in a row); see plotspec.layout_options
    width_mm, height_mm   page size; defaults per layout (plotspec.DEFAULT_SIZES_MM)
    format     'PNG', 'SVG' or 'PDF' (default: from the file suffix); dpi is used for PNG
    shared_legend  one legend for the figure (right) instead of one per panel
    title      draw the plot title above the panels
    Returns ``path``.
    """
    import matplotlib as mpl
    from matplotlib.figure import Figure

    if format is None:
        format = str(path).rsplit('.', 1)[-1]
    format = format.lower()
    default_w, default_h = DEFAULT_SIZES_MM.get(layout, DEFAULT_SIZES_MM['one'])
    width_mm, height_mm = float(width_mm or default_w), float(height_mm or default_h)
    indices, rows, cols = layout_panels(spec, layout, panel)
    rc = {'font.size': FONT_PT, 'axes.titlesize': TITLE_PT, 'axes.labelsize': FONT_PT,
          'xtick.labelsize': FONT_PT, 'ytick.labelsize': FONT_PT, 'legend.fontsize': FONT_PT,
          'figure.titlesize': TITLE_PT, 'svg.fonttype': 'none', 'pdf.fonttype': 42,
          'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6}
    with mpl.rc_context(rc):
        fig = Figure(figsize=(width_mm / MM_PER_INCH, height_mm / MM_PER_INCH), layout='constrained')
        fig.get_layout_engine().set(w_pad=1 / 72, h_pad=1 / 72, wspace=0.02, hspace=0.02)
        axes = fig.subplots(rows, cols, squeeze=False).ravel()
        handles = {}
        lw = line_width * PT_PER_PX
        for ax, index in zip(axes, indices):
            _draw_axes(fig, ax, spec.panels[index], theme, lw, handles, show_title=len(indices) > 1,
                       chars=max(12, int(width_mm / cols / 1.8)))
            if not shared_legend:
                _axes_legend(ax, spec, lw)
        for ax in axes[len(indices):]:
            ax.set_axis_off()
        if shared_legend:
            items = [handles[k][1:] for k in sorted(handles, key=lambda k: handles[k][0])] + _style_handles(spec, lw)
            if items:
                fig.legend([h for h, _ in items], [label for _, label in items], loc='outside right upper',
                           frameon=False, handlelength=2.2, borderaxespad=0.2)
        if title and spec.title:
            heading = spec.title
            if len(indices) == 1 and spec.panels[indices[0]].title:
                heading = f'{spec.title} · {spec.panels[indices[0]].title}'
            fig.suptitle(heading, wrap=True)
        if spec.message:
            fig.text(0.5, 0.5, spec.message, ha='center', va='center', color='#757575')
        fig.savefig(path, format=format, dpi=dpi)
    return path


def _draw_axes(fig, ax, panel, theme, lw, handles, show_title, chars=40):
    import textwrap

    from matplotlib.transforms import blended_transform_factory

    if show_title and panel.title:
        ax.set_title(panel.title)
    if panel.empty:
        ax.set_axis_off()
        ax.text(0.5, 0.5, textwrap.fill(panel.empty_message, chars), ha='center', va='center',
                transform=ax.transAxes, color='#757575')
        return
    if panel.heatmap is not None:
        _draw_heatmap(fig, ax, panel)
        return
    for band in panel.bands:
        ax.axvspan(band.x0, band.x1, color=BAND_COLOUR, alpha=0.2, lw=0, zorder=0)
    for s in panel.series:
        colour, dash = resolve_style(s, theme)
        marker, filled = _MPL_MARKERS.get(s.marker, (None, True)) if s.marker else (None, True)
        kwargs = {'color': colour, 'lw': lw, 'ls': _MPL_DASHES.get(dash, '-') if s.lines else 'none'}
        if marker:
            kwargs.update(marker=marker, ms=4 if not filled else 2.5, mfc=colour if filled else 'none',
                          mew=max(lw * 0.6, 0.5))
        line, = ax.plot(list(s.x), [float(v) if _finite(v) else math.nan for v in s.y],
                        label=s.name if s.legend else f'_{s.name}', **kwargs)
        key = s.key or s.name
        if s.legend and key not in handles:  # legend order: palette order, then measured data
            order = (int(s.measured), s.colour if isinstance(s.colour, int) else 10**6, len(handles))
            handles[key] = (order, line, s.name)
    xs = [float(v) for s in panel.series for v in s.x if _finite(v)]
    xs += [v.x for v in panel.vlines] + [b.x0 for b in panel.bands] + [b.x1 for b in panel.bands]
    if xs and max(xs) > min(xs):  # whole grid, also where the values are gaps
        pad = 0.02 * (max(xs) - min(xs))
        ax.set_xlim(min(xs) - pad, max(xs) + pad)
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    for vline in panel.vlines:
        ax.axvline(vline.x, ls=_MPL_DASHES.get(vline.dash, '--'), color=NEUTRAL, lw=0.8)
        if vline.label:
            ax.annotate(vline.label, (vline.x, 0.98), xycoords=trans, xytext=(-2, 0), textcoords='offset points',
                        rotation=90, ha='right', va='top', color=NEUTRAL, fontsize=FONT_PT)
    if panel.log:
        ax.set_yscale('log')
    ax.set_xlabel(panel.x_label)
    ax.set_ylabel(panel.y_label)
    ax.grid(alpha=0.25, lw=0.5)


def _draw_heatmap(fig, ax, panel):
    import numpy as np

    h = panel.heatmap
    z = np.ma.masked_invalid(np.asarray(h.z, dtype=float))
    mesh = ax.pcolormesh(np.arange(len(h.x) + 1), np.arange(len(h.y) + 1), z, cmap=h.colourscale.lower())
    ax.set_xticks(np.arange(len(h.x)) + 0.5, h.x, rotation=90 if len(h.x) > 8 else 0)
    ax.set_yticks(np.arange(len(h.y)) + 0.5, h.y)
    ax.set_xlabel(panel.x_label)
    ax.set_ylabel(panel.y_label)
    fig.colorbar(mesh, ax=ax, label=h.colourbar)


def _style_handles(spec, lw):
    from matplotlib.lines import Line2D

    return [(Line2D([], [], color=NEUTRAL, lw=lw, ls=_MPL_DASHES.get(s.dash, '-')), s.label)
            for s in spec.line_styles]


def _axes_legend(ax, spec, lw):
    if not ax.axison:
        return
    seen, items = set(), []
    for line in ax.get_lines():
        label = line.get_label()
        if label and not label.startswith('_') and label not in seen:
            seen.add(label)
            items.append((line, label))
    items += _style_handles(spec, lw)
    if items:
        ax.legend([h for h, _ in items], [label for _, label in items], frameon=False, handlelength=2.2)


def line_style_text(style):
    """'— Result A', '- - Result B' for text legends and banners."""
    return f'{DASH_SYMBOLS.get(style.dash, "—")} {style.label}'
