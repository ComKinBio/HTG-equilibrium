"""Deviation tables and plot: Julia option 2 vs pychnosz (SUPCRT92).

    uv run python julia/spikes/hkf_water/analyse.py

Reads pychnosz_reference.json and julia_results.json, prints Markdown tables
(used in README.md), writes deviations.csv and deviations.png.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
R = 8.314462618
WATER = ['rho', 'epsilon', 'QBorn', 'XBorn', 'YBorn', 'alpha', 'daldT', 'beta']
RHO_CUT = 350.0   # kg/m3, the app's aqueous validity check (phases.py, sweep.py)


def load():
    ref = json.loads((HERE / 'pychnosz_reference.json').read_text())
    jl = json.loads((HERE / 'julia_results.json').read_text())
    ions = ref['ions']
    rows = []
    for p, q in zip(ref['points'], jl['points']):
        assert (p['T_C'], p['P_bar']) == (q['T_C'], q['P_bar'])
        T = p['T_C'] + 273.15
        row = {'T_C': p['T_C'], 'P_bar': p['P_bar']}
        for k in WATER:
            s = p['water'][k]
            row[f'{k}_supcrt'] = s
            row[f'{k}_rel_iapws'] = (q['water_iapws'][k] - s) / abs(s)
            row[f'{k}_rel_check'] = (q['water_check'][k] - s) / abs(s)
        for ion in ions:
            s = p['mu0_ion'][ion]
            m = q['mu0_ion'][ion]
            row[f'{ion}_mu0'] = s
            for v in ('check', 'iapws', 'iapws_supref'):
                row[f'{ion}_dG_{v}'] = m[v] - s            # J/mol
                row[f'{ion}_rel_{v}'] = (m[v] - s) / abs(s)
                row[f'{ion}_dGRT_{v}'] = (m[v] - s) / (R * T)
        rows.append(row)
    return pd.DataFrame(rows), ions, ref, jl


def region(df):
    main = (df.T_C >= 300) & (df.T_C <= 600) & (df.P_bar >= 200) & (df.P_bar <= 300) \
        & (df.T_C % 10 == 0) & (df.P_bar % 25 == 0)
    return {
        'plan grid 300-600 C x 200-300 bar (all)': main,
        'plan grid, rho_SUPCRT >= 350 (where the app uses HKF)': main & (df.rho_supcrt >= RHO_CUT),
        'plan grid, rho_SUPCRT < 350 (app refuses electrolyte)': main & (df.rho_supcrt < RHO_CUT),
        '100-350 C at 240 bar (subcritical, app range)': (df.P_bar == 240) & (df.T_C <= 350),
        '355-390 C at 240 bar, rho >= 350 (near-critical, app range)':
            (df.P_bar == 240) & (df.T_C > 350) & (df.rho_supcrt >= RHO_CUT),
        'near-critical scan 360-420 C, 225-300 bar, rho >= 350':
            (df.T_C >= 360) & (df.T_C <= 420) & (df.rho_supcrt >= RHO_CUT),
        'near-critical scan, rho >= 350, without SUPCRT92 non-converged points':
            (df.T_C >= 360) & (df.T_C <= 420) & (df.rho_supcrt >= RHO_CUT) & ~stuck(df),
    }


def stuck(df):
    """SUPCRT92 (H2O92, LVS region) returns a fixed density of 420.04 or 199.78
    kg/m3 in narrow bands near the pseudo-critical line when its density
    iteration fails; P recomputed from that density is off by up to ~2 bar."""
    return ((df.rho_supcrt - 420.04).abs() < 0.005) | ((df.rho_supcrt - 199.78).abs() < 0.005)


def fmt(x):
    return f'{x:.1e}' if (abs(x) < 1e-2 or abs(x) >= 1e4) and x != 0 else f'{x:.3g}'


def table_water(df, mask, kind='iapws'):
    lines = ['| quantity | max abs rel dev | median abs rel dev | at max (T C, P bar) |',
             '|---|---|---|---|']
    sub = df[mask]
    for k in WATER:
        a = sub[f'{k}_rel_{kind}'].abs()
        i = a.idxmax()
        lines.append(f'| {k} | {fmt(a.max())} | {fmt(a.median())} | {sub.T_C[i]:g}, {sub.P_bar[i]:g} |')
    return '\n'.join(lines)


def table_ions(df, mask, ions, kind='iapws'):
    lines = ['| ion | max abs dG (J/mol) | median abs dG (J/mol) | max abs dG/RT | max abs rel | at max (T C, P bar) |',
             '|---|---|---|---|---|---|']
    sub = df[mask]
    for ion in ions:
        a = sub[f'{ion}_dG_{kind}'].abs()
        i = a.idxmax()
        lines.append(f'| {ion} | {fmt(a.max())} | {fmt(a.median())} | '
                     f'{fmt(sub[f"{ion}_dGRT_{kind}"].abs().max())} | {fmt(sub[f"{ion}_rel_{kind}"].abs().max())} | '
                     f'{sub.T_C[i]:g}, {sub.P_bar[i]:g} |')
    return '\n'.join(lines)


def plot(df, ions):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    for P, colour in zip((200, 225, 240, 250, 275, 300), plt.cm.viridis(np.linspace(0, 0.9, 6))):
        s = df[df.P_bar == P].sort_values('T_C')
        ok = s.rho_supcrt >= RHO_CUT
        for ax, y in zip(axes[:2], ('rho_rel_iapws', 'epsilon_rel_iapws')):
            ax.plot(s.T_C, 100 * s[y], color=colour, lw=0.8, alpha=0.4)
            ax.plot(s.T_C[ok], 100 * s[y][ok], color=colour, lw=1.8, label=f'{P} bar')
        worst = s[[f'{ion}_dG_iapws_supref' for ion in ions]].abs().max(axis=1)
        axes[2].semilogy(s.T_C, worst, color=colour, lw=0.8, alpha=0.4)
        axes[2].semilogy(s.T_C[ok], worst[ok], color=colour, lw=1.8, label=f'{P} bar')
    axes[0].set_ylabel('density, IAPWS-95 vs SUPCRT92 (%)')
    axes[1].set_ylabel('dielectric constant, JN91(IAPWS-95) vs SUPCRT92 (%)')
    axes[2].set_ylabel('max |dG°| over the 8 ions (J/mol)')
    for ax in axes:
        ax.set_xlabel('T (°C)')
        ax.axvspan(360, 420, color='0.9', zorder=0)
        ax.grid(alpha=0.3)
    axes[2].axhline(1e-10 * 3e5, color='k', ls=':', lw=1)
    axes[2].text(105, 1e-10 * 3e5 * 1.6, 'rtol 1e-10 on |mu0| ~ 3e5 J/mol', fontsize=8)
    axes[0].legend(fontsize=8)
    fig.suptitle('Option 2 vs pychnosz. Thick: water density >= 350 kg/m3 (app uses HKF); '
                 'thin: below the cut-off. Shaded: near-critical scan.', fontsize=10)
    fig.tight_layout()
    fig.savefig(HERE / 'deviations.png', dpi=110)


def main():
    df, ions, ref, jl = load()
    df.to_csv(HERE / 'deviations.csv', index=False, float_format='%.10g')
    regions = region(df)
    print('## Implementation check (SUPCRT92 density/alpha/beta fed into the Julia code)\n')
    allpts = df.T_C > 0
    print(table_water(df, allpts, 'check'))
    print()
    print(table_ions(df, allpts, ions, 'check'))
    for name, mask in regions.items():
        print(f'\n## {name}  (n = {int(mask.sum())})\n')
        print(table_water(df, mask))
        print('\nIons, option 2 with all water properties from IAPWS-95 (incl. Tr,Pr):\n')
        print(table_ions(df, mask, ions, 'iapws'))
        print('\nIons, option 2 keeping SUPCRT92 epsilon(Tr,Pr), Y(Tr,Pr):\n')
        print(table_ions(df, mask, ions, 'iapws_supref'))
    print('\n## SUPCRT92 non-converged grid points')
    print(df[stuck(df)][['T_C', 'P_bar', 'rho_supcrt', 'rho_rel_iapws']].to_string(index=False))
    print('\n## Tr,Pr water (25 C, 1 bar)')
    for k in ('rho', 'epsilon', 'YBorn'):
        print(k, ref['TrPr'][k], jl['TrPr_iapws'][k], (jl['TrPr_iapws'][k] - ref['TrPr'][k]) / ref['TrPr'][k])
    print('\n## Points at 240 bar (and thermo.json states)')
    cols = ['T_C', 'P_bar', 'rho_supcrt', 'rho_rel_iapws', 'epsilon_rel_iapws'] + \
        [f'{i}_dG_iapws_supref' for i in ('Na+', 'HCO3-', 'CO3^2-', 'Mg2+')]
    states = [(200, 240), (300, 240), (350, 250), (375, 275), (450, 240), (550, 275)]
    pick = df[((df.P_bar == 240) & df.T_C.isin([100, 150, 200, 250, 300, 340, 350, 360, 370, 375, 380, 385, 389.6, 390]))
              | df[['T_C', 'P_bar']].apply(tuple, axis=1).isin(states)]
    print(pick[cols].to_string(index=False, float_format=lambda x: f'{x:.4g}'))
    # Where does the app's 350 kg/m3 cut-off land in the two models?
    print('\n## rho = 350 kg/m3 crossing (T C) per pressure: SUPCRT92 vs IAPWS-95')
    for P in sorted(df.P_bar.unique()):
        s = df[(df.P_bar == P)].sort_values('T_C')
        rho_i = s.rho_supcrt * (1 + s.rho_rel_iapws)
        def cross(r):
            idx = np.where((r.values[:-1] >= RHO_CUT) & (r.values[1:] < RHO_CUT))[0]
            if not len(idx):
                return None
            j = idx[0]
            t0, t1, r0, r1 = s.T_C.values[j], s.T_C.values[j + 1], r.values[j], r.values[j + 1]
            return t0 + (RHO_CUT - r0) * (t1 - t0) / (r1 - r0)
        print(P, cross(s.rho_supcrt), cross(rho_i))
    plot(df, ions)


if __name__ == '__main__':
    main()
