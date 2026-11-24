"""Write test/data/formulas_python.json: htg.formulas outputs that test_formulas.jl compares with.

    uv run python julia/HTGEquilibrium/test/dump_formulas.py
"""
from pathlib import Path

import json
from htg.formulas import parse_formula, split_charge, formula_summary, hill_formula, parse_element_list, is_plain_formula, molar_mass
cases = ['CH3OH', 'CO3^2-', 'CaSO4*2H2O', 'Na+', 'CO3--', 'Mg+2', 'SO4^2-', 'Ca(OH)2', 'HCO3-', 'NH4+', 'H2O(aq)',
         'CH4(g)', 'K4[Fe(CN)6]', 'C6H12O6', 'Fe2(SO4)3', 'CaCO3·0.5H2O', 'C1.2H2.5O0.7', 'Cl-', 'PO4^3-', 'VO2+',
         'Xx', 'C2(', 'C)2', '3', 'Ca(OH', '', '  ', 'H2O^-', 'Fe^+3', 'Ca+++', '⁺Na⁺', 'Na⁺', 'CO₂', 'c2h6', 'C2H6O', 'C 2 H 6',
         'C999999.5H2', 'C0.000099999995O2']  # %g: Julia's Printf differs here
charged = [('VO2+', 1), ('VO2+', None), ('SO4-2', -2), ('NH4+', 1), ('Fe+++', 3), ('Mn^2+', 2), ('HPO4^-2', -2), ('Na', 1)]
out = {'parse': [], 'charged': [], 'summary': [], 'hill': [], 'elements': [], 'plain': []}
def run(f, *a):
    try:
        r = f(*a)
        return {'ok': r}
    except ValueError as e:
        return {'error': str(e)}
for c in cases:
    r = run(parse_formula, c)
    if 'ok' in r:
        el, z = r['ok']; r = {'elements': [[k, float(v)] for k, v in el.items()], 'charge': z, 'mass': molar_mass(el)}
    out['parse'].append([c, r])
    out['summary'].append([c, run(formula_summary, c)])
    out['plain'].append([c, is_plain_formula(c)])
for c, z in charged:
    out['charged'].append([c, z, list(split_charge(c, z)), run(formula_summary, c, z)])
for el in [{'O': 1, 'C': 1, 'H': 4}, {'Na': 1, 'Cl': 1}, {'H': 2, 'O': 1}, {'C': 1, 'O': 2}, {'C': 2.5, 'H': 1, 'S': 0.5}, {'Ca': 1, 'C': 1, 'O': 3}, {'C': 999999.5, 'H': 1e-05}]:
    out['hill'].append([[[k, v] for k, v in el.items()], hill_formula(el)])
for t in ['C, H, O', 'c;h  n', 'Na Cl, k', 'C, Xx, Yy', '', 'C,C,H']:
    out['elements'].append([t, run(parse_element_list, t)])
OUT = Path(__file__).parent / "data" / "formulas_python.json"
OUT.write_text(json.dumps(out, ensure_ascii=False, indent=None, separators=(",", ":")) + "\n")
