"""Dump the pychnosz (SUPCRT92) values the app uses, on the comparison grid.

Run from the repository root:  uv run python julia/spikes/hkf_water/dump_pychnosz.py

Writes julia/spikes/hkf_water/pychnosz_reference.json with, for every (T, P):
the eight water properties of `standard_states.supcrt_water` and the app's ion
standard Gibbs energies `standard_states.hkf_ion_g0` (J/mol, NASA reference)
for every ion in the shipped databases. The HKF records are read exactly as
the app reads them (data/supcrt92.xml through the species registry).
"""
import json
from pathlib import Path

import numpy as np

import htg
from htg import standard_states as ss
from htg.backends import chnosz

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def grid():
    points = set()
    # Plan section 4: 300-600 C x 200-300 bar, 10 C x 25 bar.
    for t in range(300, 601, 10):
        for p in range(200, 301, 25):
            points.add((float(t), float(p)))
    # The app's subcritical electrolyte range: example cases run at 240 bar,
    # 100-390 C, with aqueous cut-offs at 380-390 C.
    for t in list(range(100, 391, 5)) + [389.6]:
        points.add((float(t), 240.0))
    # Near-critical scan at the pressures the cases and reference states use.
    for p in (225.0, 240.0, 250.0, 275.0, 300.0):
        for t in np.arange(360.0, 420.01, 1.0):
            points.add((round(float(t), 3), p))
    # tests/reference/data/thermo.json states.
    for t, p in ((200, 240), (300, 240), (350, 250), (375, 275), (450, 240), (550, 275)):
        points.add((float(t), float(p)))
    return sorted(points)


def main():
    model = htg.Model.load(ROOT / 'data')
    registry = model.registry
    table = registry.ion_hkf_table
    ions = list(table['Species'])
    hkf_records = {str(r['Species']): {k: (float(r[k]) if k not in ('Species', 'SUPCRT name') else str(r[k]))
                                       for k in ('Species', 'SUPCRT name', 'Charge', 'Gf', 'Hf', 'Sr',
                                                 'a1', 'a2', 'a3', 'a4', 'c1', 'c2', 'wref')}
                   for _, r in table.iterrows()}
    ref = chnosz().water(['rho', 'epsilon', 'YBorn', 'QBorn', 'XBorn'], T=298.15, P=1.0)
    out = {'ions': ions, 'hkf_records': hkf_records,
           'TrPr': {k: float(np.asarray(v).reshape(-1)[0]) for k, v in ref.items()},
           'points': []}
    for t_c, p_bar in grid():
        T, P = t_c + 273.15, p_bar * 1e5
        row = {'T_C': t_c, 'P_bar': p_bar, 'water': ss.supcrt_water(T, P)}
        try:
            row['mu0_ion'] = dict(zip(ions, map(float, ss.hkf_ion_g0(table, ions, T, P))))
        except Exception as exc:  # noqa: BLE001 - recorded, not fatal
            row['mu0_ion'] = None
            row['error'] = f'{type(exc).__name__}: {exc}'
        out['points'].append(row)
    path = HERE / 'pychnosz_reference.json'
    path.write_text(json.dumps(out, indent=1, allow_nan=True))
    print(f'{len(out["points"])} points -> {path}')


if __name__ == '__main__':
    main()
