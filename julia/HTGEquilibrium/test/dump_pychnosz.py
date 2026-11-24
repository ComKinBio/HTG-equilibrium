"""Reference data for the Julia water/HKF tests, written by the Python app.

Run from the repository root:

    uv run python julia/HTGEquilibrium/test/dump_pychnosz.py          # HKF fixture only
    uv run python julia/HTGEquilibrium/test/dump_pychnosz.py --grid   # + wide pychnosz grid

Writes
- test/data/hkf_ions.toml (committed, ~2 kB): the HKF records of every ion in
  the shipped databases exactly as the app's registry reads them from
  data/supcrt92.xml (`registry.ion_hkf_table`). The Julia tests need them
  until registry.jl reads supcrt92.xml itself.
- test/data/pychnosz_grid.json (with --grid; generated, not committed): the
  app's `supcrt_water` and `hkf_ion_g0` on a wide T, P grid. runtests.jl
  compares against it when the file exists.
"""
import json
import sys
from pathlib import Path

import numpy as np

import htg
from htg import standard_states as ss

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
FIELDS = ('Species', 'SUPCRT name', 'Charge', 'Gf', 'Hf', 'Sr',
          'a1', 'a2', 'a3', 'a4', 'c1', 'c2', 'wref')


def toml_value(v):
    return json.dumps(v) if isinstance(v, str) else repr(float(v))


def write_fixture(table):
    lines = ['# HKF records of the ions in the shipped databases, as read by the app',
             '# registry from data/supcrt92.xml (registry.ion_hkf_table).',
             '# Written by test/dump_pychnosz.py; do not edit by hand.', '']
    for _, row in table.iterrows():
        lines.append(f'[{json.dumps(str(row["Species"]))}]')
        lines += [f'{json.dumps(k)} = {toml_value(row[k])}' for k in FIELDS]
        lines.append('')
    path = HERE / 'data' / 'hkf_ions.toml'
    path.parent.mkdir(exist_ok=True)
    path.write_text('\n'.join(lines))
    print(f'{len(table)} ions -> {path}')


def grid():
    points = set()
    for t in range(300, 601, 10):                      # plan grid
        for p in range(200, 301, 25):
            points.add((float(t), float(p)))
    for t in list(range(100, 391, 5)) + [389.6]:       # case range at 240 bar
        points.add((float(t), 240.0))
    for p in (225.0, 240.0, 250.0, 275.0, 300.0):      # near-critical scan
        for t in np.arange(360.0, 420.01, 0.5):
            points.add((round(float(t), 3), p))
    return sorted(points)


def write_grid(table):
    ions = list(table['Species'])
    out = {'ions': ions, 'points': []}
    for t_c, p_bar in grid():
        T, P = t_c + 273.15, p_bar * 1e5
        row = {'T_C': t_c, 'P_bar': p_bar, 'water': ss.supcrt_water(T, P)}
        try:
            row['mu0_ion'] = [float(x) for x in ss.hkf_ion_g0(table, ions, T, P)]
        except Exception as exc:  # noqa: BLE001 - recorded, not fatal
            row['error'] = f'{type(exc).__name__}: {exc}'
        out['points'].append(row)
    path = HERE / 'data' / 'pychnosz_grid.json'
    path.write_text(json.dumps(out))
    print(f'{len(out["points"])} points -> {path}')


def main():
    table = htg.Model.load(ROOT / 'data').registry.ion_hkf_table
    write_fixture(table)
    if '--grid' in sys.argv[1:]:
        write_grid(table)


if __name__ == '__main__':
    main()
