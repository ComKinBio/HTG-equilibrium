"""Reference cases in the settings-dict format of the frozen legacy model.

Each case reproduces a configuration from the original case notebooks
(HTG_Cases.ipynb sections 1, 2, 4, 5 and HTG_Configurable_Case.ipynb), at a few
temperatures. ``build_cases(legacy_model, data_dir)`` returns a list of dicts:

    name            unique id
    model           'base' (no species overrides) or 'acetone' (acetone -> CH3COCH3)
    atomic_weights  atomic-weight overrides applied to the model for this solve
                    (the notebooks replaced model.AW the same way)
    settings        the solve_equilibrium settings dict
    continue_from   optional name of an earlier case whose final amounts are
                    used as initial_phase_amounts (continuation)
"""
import copy

import numpy as np
import pandas as pd

GAS8 = ['water', 'hydrogen', 'carbon dioxide', 'carbon monoxide', 'methane', 'methanol',
        'ammonia', 'hydrogen sulfide']

SALT_ELEMENTS = {
    'NaHCO3': {'Na': 1, 'H': 1, 'C': 1, 'O': 3},
    'KCl': {'K': 1, 'Cl': 1},
    'MgCO3': {'Mg': 1, 'C': 1, 'O': 3},
    'CsCl': {'Cs': 1, 'Cl': 1},
    'LiCl': {'Li': 1, 'Cl': 1},
    'MgCl2': {'Mg': 1, 'Cl': 2},
}

MINERALS = ['Artinite', 'Hydromagnesite', 'Brucite', 'Magnesite', 'Periclase', 'Halite',
            'Sylvite', 'Sodium-Oxide', 'Potassium-Oxide']


def cellulose_workbook(path, AW):
    """HTG_Cases.ipynb cell 4: dry feed, atomic weights and added water from Cellulose.xlsx."""
    sheet = pd.read_excel(path, sheet_name='Input', header=None)
    feed = {str(sheet.iat[r, 0]).strip(): float(sheet.iat[r, 1]) for r in range(6, 11)}
    weights = {str(sheet.iat[r, 0]).strip(): float(sheet.iat[r, 2]) for r in range(6, 11)}
    solid_kg = float(sheet.iat[7, 6])
    water_kg = float(sheet.iat[6, 6]) - solid_kg
    carbon_kmol = solid_kg * feed['C'] / 100 / weights['C']
    water_per_mol_C = (water_kg / 18.01528) / carbon_kmol
    aw = {**AW, **weights}
    c_per_g = feed['C'] / 100 / aw['C']
    water = water_per_mol_C * c_per_g * (2 * aw['H'] + aw['O'])
    return feed, weights, water


def build_cases(legacy, data_dir):
    AW = dict(legacy.AW)
    cases = []

    def add(name, settings, model='base', atomic_weights=None, continue_from=None):
        cases.append({'name': name, 'model': model, 'atomic_weights': atomic_weights or {},
                      'settings': copy.deepcopy(settings), 'continue_from': continue_from})

    # --- HTG_Cases section 1: cellulose (molecular aqueous model, 90 % gas carbon) ---
    feed, weights, water = cellulose_workbook(data_dir / 'Cellulose.xlsx', AW)
    cellulose = {
        'pressure_bar': 275.0, 'feed_element_wt_percent': feed, 'feed_species_moles': None,
        'water_g_per_g_dry_feed': water, 'added_salts': {},
        'phase_species': {'gas': GAS8, 'aqueous': GAS8, 'organic liquid': ['benzene', 'toluene'],
                          'solid': ['carbon', 'cellulose', 'Lignin1']},
        'gas_carbon_fraction': 0.9, 'solid_carbon_fraction': None, 'solid_carbon_target': 'CHAR',
        'minimum_amounts': {('gas', 'water'): 1e-8}, 'initial_phase_amounts': None,
        'aqueous_model': 'molecular', 'aqueous_water_floor': 0.0,
        'liquid_initialization': False, 'gas_water_initialization': True,
        'optimizer_ftol': 1e-11, 'max_iterations': 1200,
    }
    for T in (400.0, 500.0):
        add(f'cellulose_{T:g}C', {**cellulose, 'temperature_C': T}, atomic_weights=weights)

    # --- HTG_Cases section 2: Yakaboylu et al. Figures 7-8 (electrolyte, 390 C cut-off) ---
    fig78 = {
        'pressure_bar': 240.0, 'aqueous_cutoff_C': 390.0,
        'feed_element_wt_percent': {e: n * AW[e] / 25000 * 100 for e, n in
                                    dict(C=632.65, H=1093.28, O=837.22, N=47.56, S=3.95).items()},
        'feed_species_moles': None,
        'water_g_per_g_dry_feed': 5550.84 * (2 * AW['H'] + AW['O']) / 25000,
        'added_salts': {},
        'phase_species': {'gas': GAS8 + ['nitrogen'], 'aqueous': list(GAS8),
                          'organic liquid': ['benzene', 'toluene'], 'solid': ['carbon']},
        'gas_carbon_fraction': None, 'solid_carbon_fraction': None, 'solid_carbon_target': 'CHAR',
        'aqueous_model': 'electrolyte', 'aqueous_water_floor': 1e-8,
        'liquid_initialization': True, 'initial_gas_water_fraction': 1e-6,
        'initial_phase_amounts': None, 'minimum_amounts': {},
        'optimizer_ftol': 1e-11, 'optimizer_eps': 1e-8, 'max_iterations': 10000,
    }
    for T in (300.0, 380.0, 450.0):
        add(f'fig78_{T:g}C', {**fig78, 'temperature_C': T})
    add('fig78_340C_continued', {**fig78, 'temperature_C': 340.0}, continue_from='fig78_300C')
    for T in (300.0, 450.0):
        add(f'fig78_{T:g}C_ipopt', {**fig78, 'temperature_C': T, 'optimizer': {'method': 'ipopt'}})
    add('fig78_450C_ipopt_multistart',
        {**fig78, 'temperature_C': 450.0, 'optimizer': {'method': 'ipopt', 'multistart': 2}})
    add('fig78_450C_slsqp_exact',
        {**fig78, 'temperature_C': 450.0, 'optimizer': {'method': 'slsqp'}})
    # Allocation constraints (solid and organic-liquid carbon) and a gas-carbon minimum.
    add('fig78_450C_allocations', {**fig78, 'temperature_C': 450.0, 'solid_carbon_fraction': 0.1,
                                   'phase_carbon_fractions': {'organic liquid': 0.05}})
    add('fig78_450C_all_solids', {**fig78, 'temperature_C': 450.0, 'solid_carbon_fraction': 0.2,
                                  'solid_carbon_target': 'all solids'})
    add('fig78_450C_gas_min', {**fig78, 'temperature_C': 450.0, 'gas_carbon_min_fraction': 0.95})
    # Dry feed given as species moles; aqueous disabled everywhere.
    add('cellulose_species_feed_450C', {**fig78, 'temperature_C': 450.0, 'aqueous_enabled': False,
                                        'feed_element_wt_percent': {},
                                        'feed_species_moles': {'cellulose': 1.0},
                                        'water_g_per_g_dry_feed': 4.0,
                                        'phase_species': {**fig78['phase_species'],
                                                          'solid': ['carbon', 'cellulose']}})

    # --- HTG_Cases section 4: full_text (ions, minerals, 380 C cut-off, exact 80 % gas C) ---
    full_text_aw = {'C': 12.011, 'H': 1.008, 'O': 15.999, 'N': 14.007, 'S': 32.06,
                    'Na': 22.989769, 'K': 39.0983, 'Cs': 132.90545196, 'Li': 6.94,
                    'Mg': 24.305, 'Cl': 35.45}
    full_text = {
        'pressure_bar': 240.0, 'aqueous_cutoff_C': 380.0, 'aqueous_enabled': True,
        'ion_product_species': {'gas': [], 'solid': list(MINERALS)},
        'feed_element_wt_percent': {'C': 30.3950366, 'H': 4.40810496, 'O': 53.57873112,
                                    'N': 2.66469168, 'S': 0.506548},
        'feed_species_moles': None, 'water_g_per_g_dry_feed': 3.999935304,
        'added_salts': {name: {'mol_per_g_dry_feed': 0.0, 'elements': SALT_ELEMENTS[name]}
                        for name in ('NaHCO3', 'KCl', 'CsCl', 'LiCl', 'MgCl2')},
        'phase_species': {'gas': GAS8 + ['nitrogen'],
                          'aqueous': GAS8 + ['Na+', 'HCO3-', 'CO3^2-', 'K+', 'Cs+', 'Cl-',
                                             'Li+', 'Mg2+'],
                          'organic liquid': ['benzene', 'toluene'], 'solid': ['carbon']},
        'gas_carbon_fraction': 0.8, 'phase_carbon_fractions': {}, 'solid_carbon_fraction': None,
        'solid_carbon_target': 'CHAR', 'aqueous_model': 'electrolyte',
        'aqueous_water_floor': 1e-6, 'minimum_amounts': {},
        'initial_gas_water_fraction': 1e-6, 'initial_phase_amounts': None,
        'optimizer_ftol': 1e-9, 'optimizer_eps': 1e-7, 'max_iterations': 5000,
    }
    liquid = {'liquid_initialization': True, 'gas_water_initialization': False}
    vapor = {'liquid_initialization': False, 'gas_water_initialization': True}
    add('full_text_250C', {**full_text, **liquid, 'temperature_C': 250.0}, atomic_weights=full_text_aw)
    add('full_text_400C', {**full_text, **vapor, 'temperature_C': 400.0}, atomic_weights=full_text_aw)
    add('full_text_250C_ipopt', {**full_text, **liquid, 'temperature_C': 250.0,
                                 'optimizer': {'method': 'ipopt'}}, atomic_weights=full_text_aw)
    # Section 5: at least 90 % gas carbon instead of exactly 80 %, with salts present.
    salty = copy.deepcopy(full_text)
    salty['added_salts']['NaHCO3']['mol_per_g_dry_feed'] = 1e-5
    salty['added_salts']['KCl']['mol_per_g_dry_feed'] = 1e-5
    add('full_text_gas_min_300C', {**salty, **liquid, 'temperature_C': 300.0,
                                   'gas_carbon_fraction': None, 'gas_carbon_min_fraction': 0.9},
        atomic_weights=full_text_aw)

    # --- HTG_Configurable_Case.ipynb defaults (acetone override, salts, organics, 385 C cut-off) ---
    organics = ['benzene', 'toluene', 'methanol', 'ethanol', 'acetone']
    neutrals = ['water', 'hydrogen', 'carbon dioxide', 'carbon monoxide', 'methane', 'ammonia',
                'hydrogen sulfide', 'nitrogen']
    configurable = {
        'pressure_bar': 240.0,
        'feed_element_wt_percent': {'C': 30.3950366, 'H': 4.40810496, 'O': 53.57873112,
                                    'N': 2.66469168, 'S': 0.506548},
        'feed_species_moles': None, 'water_g_per_g_dry_feed': 3.999935304,
        'added_salts': {name: {'mol_per_g_dry_feed': amount, 'elements': SALT_ELEMENTS[name]}
                        for name, amount in [('NaHCO3', 1e-5), ('KCl', 1e-5), ('MgCO3', 1e-5),
                                             ('CsCl', 0.0), ('LiCl', 0.0), ('MgCl2', 0.0)]},
        'phase_species': {
            'gas': neutrals + organics,
            'aqueous': neutrals + organics + ['Na+', 'K+', 'Mg2+', 'Cl-', 'HCO3-', 'CO3^2-'],
            'organic liquid': list(organics),
            'solid': ['carbon'] + MINERALS[:-2] + ['Sodium-Oxide', 'Potassium-Oxide'],
        },
        'aqueous_cutoff_C': 385.0, 'aqueous_enabled': True,
        'gas_carbon_fraction': None, 'gas_carbon_min_fraction': None,
        'solid_carbon_fraction': None, 'solid_carbon_target': 'CHAR', 'phase_carbon_fractions': {},
        'aqueous_model': 'electrolyte', 'aqueous_water_floor': 1e-6, 'minimum_amounts': {},
        'initial_phase_amounts': None, 'initial_gas_water_fraction': 1e-6,
        'max_iterations': 5000, 'optimizer_ftol': 1e-9, 'optimizer_eps': 1e-7,
    }
    add('configurable_300C', {**configurable, **liquid, 'temperature_C': 300.0}, model='acetone')
    add('configurable_450C', {**configurable, **vapor, 'temperature_C': 450.0}, model='acetone')
    add('configurable_300C_ipopt', {**configurable, **liquid, 'temperature_C': 300.0,
                                    'optimizer': {'method': 'ipopt'}}, model='acetone')

    # --- Expected failures: the new code must fail in the same way ---
    add('fail_no_carrier', {**configurable, **liquid, 'temperature_C': 300.0,
                            'added_salts': {'LiCl': {'mol_per_g_dry_feed': 1e-5,
                                                     'elements': SALT_ELEMENTS['LiCl']}},
                            'phase_species': {**configurable['phase_species'],
                                              'aqueous': neutrals + organics}}, model='acetone')
    add('fail_infeasible_fractions', {**fig78, 'temperature_C': 450.0, 'gas_carbon_fraction': 0.7,
                                      'solid_carbon_fraction': 0.5})
    return cases


# Thermodynamic states for the composition-level reference: (T [C], P [bar], electrolyte?)
STATES = [(200.0, 240.0, True), (300.0, 240.0, True), (350.0, 250.0, True),
          (375.0, 275.0, False), (450.0, 240.0, False), (550.0, 275.0, False)]

FLUID_PHASES = {
    'gas': (GAS8 + ['nitrogen'], 'vapour'),
    'organic liquid': (['benzene', 'toluene', 'methanol'], 'liquid'),
    'aqueous (molecular)': (list(GAS8), 'liquid'),
}

ELECTROLYTE_SELECTIONS = {
    'default': list(GAS8) + ['Na+', 'HCO3-', 'CO3^2-'],
    'full_text': GAS8 + ['Na+', 'HCO3-', 'CO3^2-', 'K+', 'Cs+', 'Cl-', 'Li+', 'Mg2+'],
}


def random_compositions(names, rng, count=5, ions=None, water=(2.0, 20.0)):
    """Positive amounts spanning trace to major; water first and dominant when present."""
    ions = ions or {}
    out = []
    for k in range(count):
        n = np.exp(rng.uniform(np.log(1e-9), np.log(1.0), size=len(names)))
        if k == count - 1:
            n[rng.choice(len(names), size=max(1, len(names) // 3), replace=False)] = 0.0
        if 'water' in names:
            n[names.index('water')] = rng.uniform(*water)
        if ions:
            # Make the ionic part electroneutral by adjusting the first cation/anion.
            z = np.array([ions.get(s, 0) for s in names], dtype=float)
            imbalance = float(z @ n)
            fix = next(i for i, s in enumerate(names) if z[i] * imbalance < 0)
            n[fix] += -imbalance / z[fix]
        out.append({s: float(v) for s, v in zip(names, n)})
    return out
