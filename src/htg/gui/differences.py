"""Settings that differ between runs, with the setting names of the Case tab.

Case files are flattened to dotted keys ('solver.multistart', 'feed.solids') and
each key is shown with the label the Case tab uses for it ('Initial guesses',
'Solids'). Both the old and the version-2 case format are covered; a key without
a label is shown as it is. Temperatures and pressures are shown in the chosen
display units, carbon fractions in %.
"""
import tomllib
from collections.abc import Mapping

import pandas as pd

from .. import units

PHASE_NAMES = {'gas': 'Gas', 'aqueous': 'Aqueous', 'solid': 'Solid', 'organic_liquid': 'Organic liquid'}
CARBON_BOUNDS = {'exact': 'exactly', 'min': 'at least', 'max': 'at most'}

# case key -> Case-tab label. '{T}' and '{P}' stand for the temperature and pressure unit;
# '[%]' marks carbon fractions (shown in %).
SETTING_LABELS = {
    'format': 'Case file format',
    'conditions.temperature_C': 'Temperatures [{T}]',
    'conditions.pressure_bar': 'Pressures [{P}]',
    'conditions.points': 'List of (T, P) points',
    # feed, version 2
    'feed.total_kg': 'Total mass [kg]',
    'feed.solids': 'Solids',
    'feed.ash': 'Ash',
    # feed, old format
    'feed.water_g_per_g_dry': 'Added water [g / g dry feed]',
    'feed.moisture_wt_percent': 'Moisture [% of wet feed]',
    'feed.ash_wt_percent': 'Ash [mass %, dry]',
    'feed.basis': 'Analysis basis',
    'feed.unmodelled_wt_percent': 'Unmodelled [mass %]',
    # aqueous phase
    'aqueous.enabled': 'Include an aqueous phase',
    'aqueous.model': 'Aqueous model',
    'aqueous.cutoff_C': 'Remove at and above [{T}]',
    'aqueous.water_floor': 'Smallest amount of aqueous water [mol / mol C]',
    # carbon split, old format
    'constraints.gas_carbon_fraction': 'Carbon split: gas, exactly [%]',
    'constraints.gas_carbon_min_fraction': 'Carbon split: gas, at least [%]',
    'constraints.solid_carbon_fraction': 'Carbon split: solid, exactly [%]',
    'constraints.solid_carbon_target': 'Solid carbon means',
    # starting point and solver
    'initial_guess.initialization': 'Initial composition',
    'initial_guess.initial_gas_water_fraction': 'Water in the gas at the start [fraction]',
    'solver.method': 'Method',
    'solver.multistart': 'Initial guesses',
    'solver.starts': 'Initial guesses',
    'solver.seed': 'Random seed for starting points',
    'solver.stability_test': 'Phase stability test',
    'solver.max_iterations': 'Maximum iterations',
    'solver.energy_scale': 'Energy scale [J / mol C]',
    'solver.fd_ftol': 'Finite-difference tolerance',
    'solver.fd_eps': 'Finite-difference step',
    'solver.distinct_composition_tol': 'Different results: composition [mol / mol C]',
    'solver.distinct_gibbs_tol': 'Different results: Gibbs energy [J / mol C]',
    'sweep.continuation': 'Continuation',
    'comparison.quantity': 'Measured data: its unit',
    'comparison.label': 'Measured data: label in the legend',
}

# key prefix -> label pattern for the rest of the key
PREFIX_LABELS = {
    'feed.elements_wt_percent.': '{0} [mass %]',
    'feed.species.': 'Feed: {0}',
    'feed.species_moles.': 'Feed: {0} [mol]',
    'feed.salts.': 'Added salt {0}',
    'phases.': 'Possible products: {phase0}',
    'constraints.phase_carbon_fractions.': 'Carbon split: {phase0_lower}, exactly [%]',
    'constraints.minimum_amounts.': 'Minimum amount: {phase0_lower}, {1} [mol / mol C]',
    'initial_guess.amounts.': 'Initial amount: {phase0_lower}, {1} [mol / mol C]',
    'solver.options.': 'Solver option {all}',
    'sweep.refine.': 'Extra temperatures: {0}',
    'atomic_weights.': 'Atomic weight of {0} [g / mol]',
}
TEMPERATURE_KEYS = {'conditions.temperature_C', 'aqueous.cutoff_C'}
PRESSURE_KEYS = {'conditions.pressure_bar'}
VALUE_TEXT = {True: 'yes', False: 'no'}
SKIPPED = ('name', 'description', 'comparison.file')


def _is_amount(value):
    """{value, unit} entries of the version-2 feed are shown as one setting."""
    return isinstance(value, Mapping) and set(value) <= {'value', 'unit'} and 'value' in value


def flatten(data, prefix=''):
    """{'solver': {'starts': 5}} -> {'solver.starts': 5}; {value, unit} stays one entry."""
    out = {}
    for key, value in data.items():
        name = f'{prefix}{key}'
        if isinstance(value, Mapping) and value and not _is_amount(value):
            out.update(flatten(value, name + '.'))
        else:
            out[name] = value
    return out


def setting_label(key, t_unit='C', p_unit='bar'):
    """The Case-tab name of a flattened case key."""
    if key in SETTING_LABELS:
        return SETTING_LABELS[key].replace('{T}', units.label(t_unit)).replace('{P}', units.label(p_unit))
    if key.startswith('constraints.carbon.'):  # version 2: constraints.carbon.<phase>.<exact|min|max>
        parts = key.split('.')[2:]
        if len(parts) == 2:
            phase = PHASE_NAMES.get(parts[0], parts[0]).lower()
            return f'Carbon split: {phase}, {CARBON_BOUNDS.get(parts[1], parts[1])} [%]'
    for prefix, pattern in PREFIX_LABELS.items():
        if key.startswith(prefix):
            rest = key[len(prefix):].split('.')
            phase = PHASE_NAMES.get(rest[0], rest[0])
            return pattern.format(*rest, *[''] * 2, phase0=phase, phase0_lower=phase.lower(), all=' '.join(rest))
    return key


def is_carbon_fraction(key):
    return setting_label(key).endswith('[%]') and key.startswith('constraints.')


def value_text(key, value, t_unit='C', p_unit='bar'):
    """A setting value as shown in the table, in the display units."""
    if value is None:
        return '(default)'
    if isinstance(value, bool):
        return VALUE_TEXT[value]
    if _is_amount(value):
        unit = value.get('unit', '')
        return f"{units.format_number(value['value'])} {unit}".strip()
    if key in TEMPERATURE_KEYS or key in PRESSURE_KEYS:
        unit = t_unit if key in TEMPERATURE_KEYS else p_unit
        values = value if isinstance(value, list) else [value]
        try:
            shown = [units.display(float(v), unit) for v in values]
        except (TypeError, ValueError):
            shown = [str(v) for v in values]
        return ', '.join(shown) if len(shown) <= 12 else f'{len(shown)} values'
    if key == 'conditions.points' and isinstance(value, list):
        return '; '.join(f'{units.display(T, t_unit)}, {units.display(P, p_unit)}' for T, P in value) \
            if len(value) <= 6 else f'{len(value)} points'
    if is_carbon_fraction(key) and isinstance(value, (int, float)):
        return units.format_number(100 * value)
    if isinstance(value, float):
        return units.format_number(value)
    if isinstance(value, list):
        return ', '.join(map(str, value)) if len(value) <= 12 else f'{len(value)} values'
    if isinstance(value, Mapping):
        return ', '.join(f'{k} = {v}' for k, v in value.items())
    return str(value)


def case_data(case):
    """A Case (or a case dict) as a plain dict of case-file keys."""
    if case is None:
        return None
    if isinstance(case, Mapping):
        return dict(case)
    return tomllib.loads(case.to_toml())


def setting_differences(cases, t_unit='C', p_unit='bar'):
    """Settings that differ between cases: columns 'Setting' and one per label.

    cases  {label: Case or case dict}; None entries (no case file) are left out
    """
    flat = {}
    for label, case in cases.items():
        data = case_data(case)
        if data is None:
            continue
        values = flatten(data)
        for key in SKIPPED:
            values.pop(key, None)
        flat[label] = values
    keys = sorted(set().union(*[set(f) for f in flat.values()])) if flat else []
    rows = []
    for key in keys:
        values = [flat[label].get(key) for label in flat]
        if any(v != values[0] for v in values[1:]):
            rows.append({'Setting': setting_label(key, t_unit, p_unit),
                         **{label: value_text(key, v, t_unit, p_unit) for label, v in zip(flat, values)}})
    return pd.DataFrame(rows, columns=['Setting', *flat])
