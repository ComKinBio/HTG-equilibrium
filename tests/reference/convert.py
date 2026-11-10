"""Translate reference entries (legacy settings format) into htg objects."""
import json
from pathlib import Path

import numpy as np

from htg import Case, Feed

DATA = Path(__file__).resolve().parent / 'data'
PHASE = {'gas': 'gas', 'aqueous': 'aqueous', 'organic liquid': 'organic_liquid', 'solid': 'solid'}
LEGACY_PHASE = {v: k for k, v in PHASE.items()}


def load(name):
    with open(DATA / f'{name}.json') as f:
        return json.load(f)


def pairs(triples):
    """[[legacy phase, species, value], ...] -> {(phase, species): value}."""
    return {(PHASE[p], s): float(v) for p, s, v in (triples or [])}


def case_from_legacy(entry):
    s = entry['settings']
    phases = {PHASE[p]: list(names) for p, names in s['phase_species'].items()}
    for phase, names in s.get('ion_product_species', {}).items():
        phases[PHASE[phase]] = list(dict.fromkeys(phases.get(PHASE[phase], []) + list(names)))
    flags = (bool(s.get('liquid_initialization', False)), bool(s.get('gas_water_initialization', False)))
    initialization = {(True, False): 'liquid', (False, True): 'vapor', (False, False): 'lp'}[flags]
    optimizer = s.get('optimizer')
    feed = Feed(elements_wt_percent=s.get('feed_element_wt_percent') or {},
                species_moles=s.get('feed_species_moles'),
                water_g_per_g_dry=float(s.get('water_g_per_g_dry_feed', 0.0)),
                salts=s.get('added_salts', {}))
    return Case(
        name=entry['name'], temperature_C=float(s['temperature_C']), pressure_bar=float(s['pressure_bar']),
        feed=feed, phases=phases, atomic_weights=dict(entry['atomic_weights']),
        aqueous_model=s.get('aqueous_model', 'electrolyte'),
        aqueous_enabled=s.get('aqueous_enabled', True), aqueous_cutoff_C=s.get('aqueous_cutoff_C'),
        aqueous_water_floor=float(s.get('aqueous_water_floor', 1e-8)),
        gas_carbon_fraction=s.get('gas_carbon_fraction'),
        gas_carbon_min_fraction=s.get('gas_carbon_min_fraction'),
        solid_carbon_fraction=s.get('solid_carbon_fraction'),
        solid_carbon_target={'CHAR': 'char', 'all solids': 'all_solids'}[s.get('solid_carbon_target', 'CHAR')],
        phase_carbon_fractions={PHASE[k]: v for k, v in s.get('phase_carbon_fractions', {}).items()},
        minimum_amounts=pairs(s.get('minimum_amounts')),
        initial_amounts=pairs(s['initial_phase_amounts']) if s.get('initial_phase_amounts') else None,
        initialization=initialization,
        initial_gas_water_fraction=float(s.get('initial_gas_water_fraction', 1e-6)),
        solver='slsqp-fd' if optimizer is None else optimizer['method'],
        multistart=0 if optimizer is None else int(optimizer.get('multistart', 0)),
        max_iterations=int(s.get('max_iterations', 5000)),
        fd_ftol=float(s.get('optimizer_ftol', 1e-9)), fd_eps=float(s.get('optimizer_eps', 1e-7)),
        energy_scale=s.get('optimizer_energy_scale'),
    )


def amounts_vector(result_pairs, pairs_order):
    return np.array([result_pairs.get(p, 0.0) for p in pairs_order])
