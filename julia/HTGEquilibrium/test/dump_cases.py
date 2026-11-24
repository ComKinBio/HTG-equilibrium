"""Write test/data/cases_python.json: what htg.case, htg.feed_v2 and htg.units give
for the case files of the app, so test_case.jl / test_feed.jl / test_units.jl can
compare the Julia port field by field.

    uv run python julia/HTGEquilibrium/test/dump_cases.py

Contents:
  species    elements of the species the cases use (from the app's registry,
             Model.load()), and whether each is a solid: the Julia tests use
             it in place of the registry
  files      every case file in cases/ (recursively), the old (format 1)
             files in tests/data/legacy_cases and the template: the loaded
             Case, its grid, chains, feed numbers, carbon-bound messages,
             feed summary, the text to_toml writes, legacy_to_v2 and the case
             converted from it
  variants   the template (or an old file) with text replaced: one setting
             or one mistake each; the same dump, or the error
  feeds      feed_v2 functions called directly on [feed] tables
  units, reprs, rounds, plus_one, arange: small helpers

Values are written as JSON with dict order kept; dicts whose keys are not
strings become {"__pairs__": [[key, value], ...]}, non-finite floats
{"__float__": "inf"}. Paths are relative to the repository root.
"""
import copy
import json
import math
import os
import tomllib
from dataclasses import fields
from decimal import Decimal  # noqa: F401  (used by feed_v2._plus_one)
from pathlib import Path

import numpy as np

import htg
from htg import units
from htg.case import (Case, Feed, Salt, carbon_bound_messages, case_from_dict, effective_carbon_bounds,
                      fluid_phase_count, legacy_to_v2, refined_temperatures)
from htg.feed_v2 import (_plus_one, ash_percent, elemental_masses, feed_fractions, feed_summary,
                         feed_v2_to_model, legacy_feed_to_v2, molar_mass_function, registry_elements)
from htg.model import Model

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).parent / 'data' / 'cases_python.json'
TEMPLATE = 'src/htg/templates/case_template.toml'
LEGACY_TEMPLATE = 'tests/data/legacy_cases/case_template.toml'


def enc(x):
    if x is None or isinstance(x, (bool, int, str)):
        return x
    if isinstance(x, float):
        return x if math.isfinite(x) else {'__float__': repr(x)}
    if isinstance(x, (np.floating, np.integer)):
        return enc(x.item())
    if isinstance(x, (list, tuple)):
        return [enc(v) for v in x]
    if isinstance(x, dict):
        if all(isinstance(k, str) for k in x):
            return {k: enc(v) for k, v in x.items()}
        return {'__pairs__': [[enc(k), enc(v)] for k, v in x.items()]}
    if isinstance(x, (Case, Feed, Salt)):
        return {f.name: enc(getattr(x, f.name)) for f in fields(x)}
    raise TypeError(f'cannot encode {x!r}')


def error(exc):
    return {'type': type(exc).__name__, 'message': str(exc)}


def attempt(f, *args, **kwargs):
    try:
        return {'ok': f(*args, **kwargs)}
    except (ValueError, KeyError, TypeError) as exc:
        return {'error': error(exc)}


# ---------------------------------------------------------------- species data
MODEL = Model.load()
ELEMENTS_OF = registry_elements(MODEL.registry)
SOLIDS = MODEL.registry.solid
FEED_SPECIES = ['water', 'cellulose', 'methanol', 'Lignin1', 'carbon']


def species_elements(name):
    return ELEMENTS_OF(name)


def molar_mass(case=None):
    return molar_mass_function(species_elements, case.atomic_weights if case else None)


def carbon_solids(case):
    return {name: bool(SOLIDS[name].get('elements', {}).get('C', 0))
            for name in case.phases.get('solid', []) if name in SOLIDS}


# ---------------------------------------------------------------- one case
def case_dump(case):
    fields_ = enc(case)
    comparison = fields_.get('comparison')
    if comparison and os.path.isabs(comparison.get('file', '')):
        comparison['file'] = os.path.relpath(comparison['file'], ROOT)
    out = {'case': fields_}
    for name in ('temperatures', 'pressures', 'grid', 'chains', 'is_single_point'):
        out[name] = enc(attempt(getattr(case, name)))
    pressures = out['pressures'].get('ok')
    temperatures = out['temperatures'].get('ok')
    if temperatures and pressures:
        point = case.at(temperatures[0], pressures[0])
        out['at'] = enc({k: getattr(point, k) for k in ('temperature_C', 'pressure_bar', 'points', 'refine')})
    solids = carbon_solids(case)
    out['carbon_solids'] = enc(solids)
    out['carbon_bound_messages'] = [carbon_bound_messages(case), carbon_bound_messages(case, solids)]
    out['effective_carbon_bounds'] = enc([effective_carbon_bounds(case.carbon_bounds, True),
                                          effective_carbon_bounds(case.carbon_bounds, False)])
    out['fluid_phase_count'] = fluid_phase_count(case.phases, case.aqueous_enabled, case.aqueous_cutoff_C,
                                                 temperatures)
    resolved = attempt(case.feed.resolved, molar_mass(case))
    numbers = {'resolved': enc(resolved)}
    if 'ok' in resolved:
        feed = resolved['ok']
        numbers['water'] = enc(attempt(feed.total_water_g_per_g_dry))
        if feed.species_moles is None:
            numbers['dry'] = enc(attempt(feed.dry_element_wt_percent, None))
    numbers['fractions'] = enc(attempt(feed_fractions, case.feed, molar_mass(case)))
    numbers['fractions_without_molar_mass'] = enc(attempt(feed_fractions, case.feed))
    out['feed'] = numbers
    if case.feed_v2 is not None:
        out['feed_summary'] = enc(attempt(feed_summary, case.feed_v2, case.atomic_weights, species_elements))
    text = attempt(case.to_toml)
    if 'ok' in text:
        text['ok'] = text['ok'].replace(str(ROOT), '<repo>')
    out['to_toml'] = text
    return out


def text_dump(text, name):
    """Everything for one case-file text."""
    data = tomllib.loads(text)
    out = {}
    try:
        case = case_from_dict(data, name=name)
    except (ValueError, KeyError, TypeError) as exc:
        out['error'] = error(exc)
    else:
        out.update(case_dump(case))
    converted = attempt(legacy_to_v2, data, molar_mass())
    out['legacy_to_v2'] = enc(converted)
    if 'ok' in converted and converted['ok'] != data:
        try:
            out['converted'] = case_dump(case_from_dict(converted['ok'], name=name))
        except (ValueError, KeyError, TypeError) as exc:
            out['converted'] = {'error': error(exc)}
    return out


def file_dump(path):
    rel = str(path.relative_to(ROOT))
    out = {'path': rel}
    try:
        case = htg.load_case(path)
    except (ValueError, KeyError, TypeError) as exc:
        out['error'] = error(exc)
        return out
    out.update(case_dump(case))
    data = tomllib.loads(path.read_text(encoding='utf-8'))
    converted = attempt(legacy_to_v2, data, molar_mass())
    out['legacy_to_v2'] = enc(converted)
    if 'ok' in converted and converted['ok'] != data:
        try:
            out['converted'] = case_dump(case_from_dict(converted['ok'], name=path.stem))
        except (ValueError, KeyError, TypeError) as exc:
            out['converted'] = {'error': error(exc)}
    return out


# ---------------------------------------------------------------- variants of the template
CARBON = ('# [constraints.carbon]\n# gas = { min = 0.9 }\n# solid = { max = 0.05 }\n'
          '# aqueous = { exact = 0.03 }         # dropped at temperatures without an aqueous phase')
GRID = ('temperature_C = [300, 350, 400, 450]   # one value, or a list for a sweep\n'
        'pressure_bar = 240.0                   # one value, or a list: every temperature at every pressure')
TEMPS = 'temperature_C = [300, 350, 400, 450]'
ELEMENTAL_FEED = ('total_kg = 5.0                              # solids + water\n'
                  'solids = { value = 20.0, unit = "%" }       # "%" of the total mass, or "kg"\n'
                  'ash = { value = 8.44, unit = "%" }          # "%" of the solids, or "kg"; leave out for no ash\n')
ELEMENTS = 'elements_wt_percent = { C = 30.4, H = 4.41, O = 53.58, N = 2.66, S = 0.51 }'
SPECIES_FEED = '[feed.species]\nwater = { value = 4.5, unit = "kg" }\ncellulose = { value = 0.5, unit = "kg" }\n'
SOLVER = 'method = "ipopt"          #'
LEGACY_WATER = 'water_g_per_g_dry = 4.0'


def carbon(table):
    return [(CARBON, '[constraints.carbon]\n' + table)]


VARIANTS = [
    ('points', TEMPLATE, [(GRID, 'points = [[300, 240], [350, 250.5], [300, 240], [400, 260]]')]),
    ('points, no continuation', TEMPLATE, [(GRID, 'points = [[300, 240], [350, 250]]'),
                                           ('continuation = true', 'continuation = false')]),
    ('pressure list', TEMPLATE, [('pressure_bar = 240.0 ', 'pressure_bar = [250, 240.0, 250.0, 225.123456789] ')]),
    ('no continuation', TEMPLATE, [('continuation = true', 'continuation = false')]),
    ('single point', TEMPLATE, [(TEMPS, 'temperature_C = 375.5')]),
    ('temperatures as text', TEMPLATE, [(TEMPS, 'temperature_C = ["300", "3_50", " 400.25 ", 450]')]),
    ('bad temperature', TEMPLATE, [(TEMPS, 'temperature_C = [300, "hot"]')]),
    ('temperature below absolute zero', TEMPLATE, [(TEMPS, 'temperature_C = [-300, 300]')]),
    ('no temperatures', TEMPLATE, [(TEMPS, 'temperature_C = []')]),
    ('bad pressure', TEMPLATE, [('pressure_bar = 240.0 ', 'pressure_bar = 0 ')]),
    ('missing pressure', TEMPLATE, [('pressure_bar = 240.0 ', '# ')]),
    ('points and temperatures', TEMPLATE, [(TEMPS, TEMPS + '\npoints = [[300, 240]]')]),
    ('bad points entry', TEMPLATE, [(GRID, 'points = [[300, 240], [350]]')]),
    ('empty points', TEMPLATE, [(GRID, 'points = []')]),
    ('refine', TEMPLATE, [('[sweep]\n', '[sweep]\nrefine = { start_C = 310, stop_C = 390, step_C = 7.5 }\n')]),
    ('refine fine', TEMPLATE, [('[sweep]\n', '[sweep]\nrefine = { start_C = 300.05, stop_C = 301, step_C = 0.1 }\n')]),
    ('refine with points', TEMPLATE, [(GRID, 'points = [[300, 240]]'),
                                      ('[sweep]\n', '[sweep]\nrefine = { start_C = 300, stop_C = 400, step_C = 5 }\n')]),
    ('cut-off inside the grid', TEMPLATE, [('cutoff_C = 390.0', 'cutoff_C = 375.5')]),
    ('aqueous disabled', TEMPLATE, [('# enabled = false', 'enabled = false')]),
    ('bad aqueous model', TEMPLATE, [('model = "electrolyte"', 'model = "ideal"')]),
    ('carbon bounds', TEMPLATE, carbon('gas = { min = 0.9 }\nsolid = { max = 0.05 }\n')),
    ('carbon bounds ints', TEMPLATE, carbon('solid = { max = 0 }\ngas = { min = 1 }\n')),
    ('carbon bound min > max', TEMPLATE, carbon('solid = { min = 0.2, max = 0.05 }\n')),
    ('carbon bound exact and min', TEMPLATE, carbon('gas = { exact = 0.5, min = 0.2 }\n')),
    ('carbon bound above 1', TEMPLATE, carbon('aqueous = { max = 1.5 }\n')),
    ('carbon bound inf', TEMPLATE, carbon('aqueous = { max = inf }\n')),
    ('carbon bound text', TEMPLATE, carbon('aqueous = { max = "half" }\n')),
    ('carbon bound bad kind', TEMPLATE, carbon('organic_liquid = { at_least = 0.1 }\n')),
    ('carbon bound unknown phase', TEMPLATE, carbon('liquid = { exact = 0.1 }\n')),
    ('aqueous bound, aqueous disabled', TEMPLATE, carbon('aqueous = { exact = 0.03 }\n')
     + [('# enabled = false', 'enabled = false')]),
    ('every phase fixed, sum 0.9', TEMPLATE, carbon('gas = { exact = 0.5 }\naqueous = { exact = 0.2 }\n'
                                                    'organic_liquid = { exact = 0.1 }\nsolid = { exact = 0.1 }\n')),
    ('every phase fixed, all solids', TEMPLATE, carbon('gas = { exact = 0.5 }\naqueous = { exact = 0.2 }\n'
                                                       'organic_liquid = { exact = 0.1 }\nsolid = { exact = 0.1 }\n')
     + [('# solid_carbon_target = "char"', 'solid_carbon_target = "all_solids"')]),
    ('every phase fixed, sum 1', TEMPLATE, carbon('gas = { exact = 0.6 }\naqueous = { exact = 0.2 }\n'
                                                  'organic_liquid = { exact = 0.1 }\nsolid = { exact = 0.1 }\n')),
    ('carbon shares above 100 %', TEMPLATE, carbon('gas = { min = 0.9 }\nsolid = { exact = 0.2 }\n')),
    ('bad solid carbon target', TEMPLATE, [('# solid_carbon_target = "char"', 'solid_carbon_target = "coke"')]),
    ('legacy carbon keys and multistart', TEMPLATE,
     [('[constraints]\n', '[constraints]\ngas_carbon_min_fraction = 0.9\nsolid_carbon_fraction = 0.05\n'
                          'phase_carbon_fractions = { organic_liquid = 0.01 }\n'),
      (SOLVER, 'multistart = 4\n' + SOLVER)]),
    ('legacy all solids', TEMPLATE, [('[constraints]\n', '[constraints]\ngas_carbon_fraction = 0.8\n'
                                                         'phase_carbon_fractions = { solid = 0.1, gas = 0.8 }\n')]),
    ('legacy solid conflict', TEMPLATE, [('[constraints]\n', '[constraints]\nsolid_carbon_fraction = 0.05\n'
                                                             'phase_carbon_fractions = { solid = 0.1 }\n')]),
    ('legacy gas conflict', TEMPLATE, [('[constraints]\n', '[constraints]\ngas_carbon_fraction = 0.8\n'
                                                           'phase_carbon_fractions = { gas = 0.7 }\n')]),
    ('legacy gas both', TEMPLATE, [('[constraints]\n', '[constraints]\ngas_carbon_fraction = 0.8\n'
                                                       'gas_carbon_min_fraction = 0.7\n')]),
    ('legacy and new gas bound', TEMPLATE, carbon('gas = { min = 0.9 }\n')
     + [('[constraints]\n', '[constraints]\ngas_carbon_fraction = 0.8\n')]),
    ('legacy and same new gas bound', TEMPLATE, carbon('gas = { exact = 0.8 }\n')
     + [('[constraints]\n', '[constraints]\ngas_carbon_fraction = 0.8\n')]),
    ('multistart, aqueous cut off everywhere', TEMPLATE, [('cutoff_C = 390.0', 'cutoff_C = 250'),
                                                          (SOLVER, 'multistart = 2\n' + SOLVER)]),
    ('starts and multistart', TEMPLATE, [(SOLVER, 'multistart = 2\nstarts = 3\n' + SOLVER)]),
    ('negative multistart', TEMPLATE, [(SOLVER, 'multistart = -1\n' + SOLVER)]),
    ('starts 3', TEMPLATE, [('# starts = 5', 'starts = 3.0')]),
    ('starts 2.5', TEMPLATE, [('# starts = 5', 'starts = 2.5')]),
    ('starts 0', TEMPLATE, [('# starts = 5', 'starts = 0')]),
    ('starts text', TEMPLATE, [('# starts = 5', 'starts = "x"')]),
    ('starts bool', TEMPLATE, [('# starts = 5', 'starts = true')]),
    ('negative tolerance', TEMPLATE, [('# distinct_gibbs_tol = 1.0', 'distinct_gibbs_tol = -1.0')]),
    ('solver options', TEMPLATE, [(SOLVER, 'options = { tol = 1e-10, print_level = 0 }\nseed = 7\n'
                                           'stability_test = false\nenergy_scale = 1e5\n' + SOLVER)]),
    ('bad solver', TEMPLATE, [('method = "ipopt"', 'method = "newton"')]),
    ('bad initialization', TEMPLATE, [('initialization = "auto"', 'initialization = "random"')]),
    ('unknown solver key', TEMPLATE, [(SOLVER, 'foo = 1\n' + SOLVER)]),
    ('unknown section', TEMPLATE, [('[sweep]\n', '[extra]\nx = 1\n\n[sweep]\n')]),
    ('unknown top-level key', TEMPLATE, [('format = 2\n', 'format = 2\nauthor = "me"\n')]),
    ('unknown phase', TEMPLATE, [('solid = ["carbon"]', 'solid = ["carbon"]\nliquid = ["water"]')]),
    ('format 3', TEMPLATE, [('format = 2\n', 'format = 3\n')]),
    ('format 1', TEMPLATE, [('format = 2\n', 'format = 1\n')]),
    ('minimum amounts', TEMPLATE, [('# [constraints.minimum_amounts]      # lower bounds, mol/mol C\n# gas = { water = 1e-8 }',
                                    '[constraints.minimum_amounts]\ngas = { water = 1e-8, methane = 2 }\n'
                                    'aqueous = { water = 0.5 }')]),
    ('initial amounts', TEMPLATE, [('\n[sweep]\n', '\n[initial_guess.amounts]\ngas = { water = 0.5, methane = 0.1 }\n\n[sweep]\n')]),
    ('comparison', TEMPLATE, [('# [comparison]\n# file = "comparison/my_data.csv"',
                               '[comparison]\nfile = "comparison/my_data.csv"\nlabel = "Paper"')]),
    ('comparison without file', TEMPLATE, [('# [comparison]\n', '[comparison]\nlabel = "Paper"\n')]),
    ('atomic weights', TEMPLATE, [('# atomic_weights = { C = 12.011, H = 1.008 }',
                                   'atomic_weights = { C = 12, H = 1.00794 }')]),
    ('salts', TEMPLATE, [('# [feed.salts.NaHCO3]\n# mol_per_g_dry_feed = 1e-5\n# elements = { Na = 1, H = 1, C = 1, O = 3 }',
                          '[feed.salts.NaHCO3]\nmol_per_g_dry_feed = 1e-5\nelements = { Na = 1, H = 1, C = 1, O = 3 }\n\n'
                          '[feed.salts.KCl]\nmol_per_g_dry_feed = 2e-5\nelements = { K = 1, Cl = 1 }')]),
    ('species feed', TEMPLATE, [(ELEMENTAL_FEED, ''), (ELEMENTS, SPECIES_FEED)]),
    ('species feed in mol', TEMPLATE, [(ELEMENTAL_FEED, ''), (ELEMENTS, '[feed.species]\nwater = { value = 250, unit = "mol" }\n'
                                                                       'cellulose = { value = 0.5, unit = "kg" }\n'
                                                                       'methanol = { value = 2, unit = "mol" }\n'
                                                                       'Lignin1 = { value = 0.125 }\n')]),
    ('species feed, unknown species', TEMPLATE, [(ELEMENTAL_FEED, ''), (ELEMENTS, SPECIES_FEED + 'unobtainium = { value = 1, unit = "kg" }\n')]),
    ('species feed, water only', TEMPLATE, [(ELEMENTAL_FEED, ''), (ELEMENTS, '[feed.species]\nwater = { value = 4.5, unit = "kg" }\n')]),
    ('species feed, no water', TEMPLATE, [(ELEMENTAL_FEED, ''), (ELEMENTS, '[feed.species]\ncellulose = { value = 4.5, unit = "kg" }\n')]),
    ('species feed, bad unit', TEMPLATE, [(ELEMENTAL_FEED, ''), (ELEMENTS, SPECIES_FEED + 'methanol = { value = 1, unit = "g" }\n')]),
    ('species feed, bad key', TEMPLATE, [(ELEMENTAL_FEED, ''), (ELEMENTS, SPECIES_FEED + 'methanol = { value = 1, units = "kg" }\n')]),
    ('species feed, negative', TEMPLATE, [(ELEMENTAL_FEED, ''), (ELEMENTS, SPECIES_FEED + 'methanol = { value = -1, unit = "kg" }\n')]),
    ('species feed, not a table', TEMPLATE, [(ELEMENTAL_FEED, ''), (ELEMENTS, SPECIES_FEED + 'methanol = 1.0\n')]),
    ('species feed and total', TEMPLATE, [('ash = { value = 8.44, unit = "%" }', ''), (ELEMENTS, SPECIES_FEED)]),
    ('mixed feed keys', TEMPLATE, [('total_kg = 5.0', 'total_kg = 5.0\nwater_g_per_g_dry = 4.0')]),
    ('elements only, format 2', TEMPLATE, [(ELEMENTAL_FEED, '')]),
    ('elements only, no format', TEMPLATE, [(ELEMENTAL_FEED, ''), ('format = 2\n', '')]),
    ('no feed', TEMPLATE, [('[feed]\n', '[feedx]\n')]),
    ('ash in kg', TEMPLATE, [('ash = { value = 8.44, unit = "%" }', 'ash = { value = 0.0844, unit = "kg" }')]),
    ('no ash', TEMPLATE, [('ash = { value = 8.44, unit = "%" }', '')]),
    ('solids in kg', TEMPLATE, [('solids = { value = 20.0, unit = "%" }', 'solids = { value = 1, unit = "kg" }')]),
    ('solids unit g', TEMPLATE, [('solids = { value = 20.0, unit = "%" }', 'solids = { value = 20.0, unit = "g" }')]),
    ('solids above total', TEMPLATE, [('solids = { value = 20.0, unit = "%" }', 'solids = { value = 6, unit = "kg" }')]),
    ('solids as a number', TEMPLATE, [('solids = { value = 20.0, unit = "%" }', 'solids = 20.0')]),
    ('solids missing', TEMPLATE, [('solids = { value = 20.0, unit = "%" }', '')]),
    ('total missing', TEMPLATE, [('total_kg = 5.0', '')]),
    ('total zero', TEMPLATE, [('total_kg = 5.0', 'total_kg = 0')]),
    ('total text', TEMPLATE, [('total_kg = 5.0', 'total_kg = "5 kg"')]),
    ('total bool', TEMPLATE, [('total_kg = 5.0', 'total_kg = true')]),
    ('total nan', TEMPLATE, [('total_kg = 5.0', 'total_kg = nan')]),
    ('ash above solids', TEMPLATE, [('ash = { value = 8.44, unit = "%" }', 'ash = { value = 2.0, unit = "kg" }')]),
    ('negative ash', TEMPLATE, [('ash = { value = 8.44, unit = "%" }', 'ash = { value = -1 }')]),
    ('elements missing', TEMPLATE, [(ELEMENTS, '')]),
    ('elements empty', TEMPLATE, [(ELEMENTS, 'elements_wt_percent = {}')]),
    ('elements sum 100', TEMPLATE, [(ELEMENTS, 'elements_wt_percent = { C = 50, H = 6.1, O = 43.9 }')]),
    ('elements as text', TEMPLATE, [(ELEMENTS, 'elements_wt_percent = { C = "30.4", H = 4.41, O = 53.58 }')]),
    ('element negative', TEMPLATE, [(ELEMENTS, 'elements_wt_percent = { C = 30.4, H = -4.41, O = 53.58 }')]),
    ('element bool', TEMPLATE, [(ELEMENTS, 'elements_wt_percent = { C = true, H = 4.41 }')]),
    ('elements zero', TEMPLATE, [(ELEMENTS, 'elements_wt_percent = { C = 0, H = 0.0 }')]),
    ('unknown element', TEMPLATE, [(ELEMENTS, 'elements_wt_percent = { C = 30.4, Xx = 4.41, O = 53.58 }')]),
    ('unknown feed key', TEMPLATE, [('total_kg = 5.0', 'total_kg = 5.0\nmass = 3')]),
    # old (format 1) files
    ('legacy unmodelled', LEGACY_TEMPLATE, [(LEGACY_WATER, LEGACY_WATER + '\nunmodelled_wt_percent = 3.0\nash_wt_percent = 5.44')]),
    ('legacy ash', LEGACY_TEMPLATE, [(LEGACY_WATER, LEGACY_WATER + '\nash_wt_percent = 8.62')]),
    ('legacy ash off', LEGACY_TEMPLATE, [(LEGACY_WATER, LEGACY_WATER + '\nash_wt_percent = 8.6')]),
    ('legacy moisture', LEGACY_TEMPLATE, [(LEGACY_WATER, LEGACY_WATER + '\nmoisture_wt_percent = 20')]),
    ('legacy moisture 100', LEGACY_TEMPLATE, [(LEGACY_WATER, LEGACY_WATER + '\nmoisture_wt_percent = 100')]),
    ('legacy dry ash free', LEGACY_TEMPLATE, [('C = 30.4, H = 4.41, O = 53.58, N = 2.66, S = 0.51', 'C = 50, H = 6, O = 44'),
                                              (LEGACY_WATER, LEGACY_WATER + '\nbasis = "dry_ash_free"\nash_wt_percent = 5.0')]),
    ('legacy dry ash free, bad sum', LEGACY_TEMPLATE, [(LEGACY_WATER, LEGACY_WATER + '\nbasis = "dry_ash_free"')]),
    ('legacy bad basis', LEGACY_TEMPLATE, [(LEGACY_WATER, LEGACY_WATER + '\nbasis = "wet"')]),
    ('legacy species moles', LEGACY_TEMPLATE, [(ELEMENTS, 'species_moles = { cellulose = 2.0, methanol = 1 }')]),
    ('legacy species moles with ash', LEGACY_TEMPLATE, [(ELEMENTS, 'species_moles = { cellulose = 2.0 }\nash_wt_percent = 3')]),
    ('legacy unknown feed key', LEGACY_TEMPLATE, [(LEGACY_WATER, LEGACY_WATER + '\nmass = 1')]),
    ('legacy refine', LEGACY_TEMPLATE, [('continuation = true', 'continuation = true\nrefine = { start_C = 340, stop_C = 370, step_C = 10 }')]),
    ('legacy refine only', LEGACY_TEMPLATE, [('continuation = true', 'refine = { start_C = 340, stop_C = 370, step_C = 10 }')]),
    ('legacy refine with points', LEGACY_TEMPLATE, [(GRID, 'points = [[300, 240], [350, 250]]'),
                                                    ('continuation = true', 'refine = { start_C = 340, stop_C = 370, step_C = 10 }')]),
    ('legacy multistart with points', LEGACY_TEMPLATE, [(GRID, 'points = [[400, 240], [450, 250]]'),
                                                        ('method = "ipopt"', 'method = "ipopt"\nmultistart = 3')]),
    ('legacy multistart', LEGACY_TEMPLATE, [('method = "ipopt"', 'method = "ipopt"\nmultistart = 3'),
                                            ('continuation = true', 'refine = { start_C = 340, stop_C = 370, step_C = 10 }')]),
]


def variant_text(base, replacements):
    text = (ROOT / base).read_text(encoding='utf-8')
    for old, new in replacements:
        assert old in text, (base, old)
        text = text.replace(old, new)
    return text


# ---------------------------------------------------------------- direct feed calls
FEED = {'total_kg': 10.0, 'solids': {'value': 15.0, 'unit': '%'}, 'ash': {'value': 8.45, 'unit': '%'},
        'elements_wt_percent': {'C': 30.4, 'H': 4.41, 'O': 53.58, 'N': 2.66, 'S': 0.51}}
SPECIES = {'species': {'water': {'value': 4.5, 'unit': 'kg'}, 'cellulose': {'value': 0.5, 'unit': 'kg'},
                       'methanol': {'value': 2.0, 'unit': 'mol'}}}
FEED_INPUTS = [
    FEED,
    {**FEED, 'solids': {'value': 1.5, 'unit': 'kg'}, 'ash': {'value': 0.12675, 'unit': 'kg'}},
    {**FEED, 'ash': {'value': 0.12675, 'unit': 'kg'}},
    {k: v for k, v in FEED.items() if k != 'ash'},
    {**FEED, 'ash': None},
    {**FEED, 'salts': {'KCl': {'mol_per_g_dry_feed': 1e-5, 'elements': {'K': 1, 'Cl': 1}}}},
    {**FEED, 'ash': {'value': 2.0, 'unit': 'kg'}},
    {**FEED, 'elements_wt_percent': {}},
    {**FEED, 'solids': {'value': 15.0, 'unit': 'g'}},
    {**FEED, 'total_kg': 0},
    {**FEED, 'solids': {'value': 12.0, 'unit': 'kg'}},
    {**FEED, 'solids': {'value': 100.0, 'unit': '%'}},
    {**FEED, 'solids': {'value': 15.0}},
    {**FEED, 'solids': {'value': 15.0, 'unit': '%', 'basis': 'wet'}},
    {**FEED, 'elements_wt_percent': {'C': 30.4, 'Xx': 1.0}},
    {**FEED, 'total_kg': '10'},
    {**FEED, 'total_kg': [10]},
    SPECIES,
    {'species': {**SPECIES['species'], 'water': {'value': 0.0, 'unit': 'kg'}}},
    {**SPECIES, 'total_kg': 5.0},
    {'species': {'cellulose': {'value': 0.5, 'unit': 'kg'}}},
]
LEGACY_FEEDS = [
    {'elements_wt_percent': {'C': 45.0, 'H': 5.0, 'O': 40.0}, 'ash_wt_percent': 10.0, 'water_g_per_g_dry': 3.0,
     'moisture_wt_percent': 20.0},
    {'elements_wt_percent': {'C': 50.0, 'H': 6.0, 'O': 44.0}, 'basis': 'dry_ash_free', 'ash_wt_percent': 5.0},
    {'elements_wt_percent': {'C': 50.0, 'H': 6.0, 'O': 44.0}},
    {'elements_wt_percent': {'C': 50, 'H': 6, 'O': 40}, 'water_g_per_g_dry': 0.1},
    {'elements_wt_percent': {'C': 30.3950366, 'H': 4.40810496, 'O': 53.57873112, 'N': 2.66469168, 'S': 0.506548},
     'water_g_per_g_dry': 3.999935304, 'ash_wt_percent': 8.44688764},
    {'elements_wt_percent': {'C': 44.2, 'H': 6.3, 'O': 49.5, 'N': 0.0, 'S': 0.0}, 'water_g_per_g_dry': 8.999660288377422},
    {'elements_wt_percent': {'C': 50.0, 'H': 6.0, 'O': 34.0}, 'ash_wt_percent': 5.0, 'unmodelled_wt_percent': 5.0},
    {'species_moles': {'cellulose': 1.0}, 'water_g_per_g_dry': 4.0},
    {'species_moles': {'cellulose': 2.0, 'methanol': 1}, 'water_g_per_g_dry': 4.0},
    {'elements_wt_percent': {'C': 50.0, 'H': 6.0}, 'water_g_per_g_dry': 1.0,
     'salts': {'KCl': {'mol_per_g_dry_feed': 1e-5, 'elements': {'K': 1, 'Cl': 1}}}},
]


def feed_calls():
    out = []
    mm = molar_mass()
    for feed in FEED_INPUTS:
        out.append({'feed': enc(feed),
                     'feed_summary': enc(attempt(feed_summary, copy.deepcopy(feed), None, species_elements)),
                     'feed_summary_no_species': enc(attempt(feed_summary, copy.deepcopy(feed))),
                     'feed_v2_to_model': enc(attempt(feed_v2_to_model, copy.deepcopy(feed))),
                     'feed_v2_to_model_mm': enc(attempt(feed_v2_to_model, copy.deepcopy(feed), molar_mass=mm)),
                     'elemental_masses': enc(attempt(elemental_masses, copy.deepcopy(feed))),
                     'ash_percent': enc(attempt(ash_percent, copy.deepcopy(feed)))})
        model = out[-1]['feed_v2_to_model']
        if 'ok' in model:
            out[-1]['feed_fractions'] = enc(attempt(feed_fractions, feed_v2_to_model(copy.deepcopy(feed)), mm))
    legacy = []
    for feed in LEGACY_FEEDS:
        legacy.append({'feed': enc(feed), 'legacy_feed_to_v2': enc(attempt(legacy_feed_to_v2, copy.deepcopy(feed))),
                       'legacy_feed_to_v2_mm': enc(attempt(legacy_feed_to_v2, copy.deepcopy(feed), mm))})
    return out, legacy


# ---------------------------------------------------------------- helpers
NUMBERS = [0.0, -0.0, 1.0, 300.0, 0.1 + 0.2, 573.15, 299.99999999999994, 1234.5678, 1234567, 1.5e-5, -40,
           2.4e7, 101325, 1e5, 24, 236.86158401184308, 123456.5, 999999.5, 9999995, 0.000123456789, 1e21,
           1e-7, 2.5e-300, 1.7976931348623157e308, 5e-324, 12345678901234567890.0, 0.5, 2.675, 1e16, 1e15,
           123456789.123, 0.0001, 0.00001, -1e-5, 3.999935304, 8.44688764, math.pi, math.e, float('inf'),
           float('-inf'), float('nan')]


def helpers():
    finite = [x for x in NUMBERS if math.isfinite(x)]
    out = {'reprs': [[enc(float(x)), repr(float(x))] for x in NUMBERS]}
    out['format_number'] = [[enc(float(x)), unit, digits, units.format_number(x, unit, digits)]
                            for x in NUMBERS for unit in (None, 'C', 'Pa', 'MPa') for digits in (6, 3)]
    out['conversions'] = [[float(x), unit, enc(units.to_file(float(x), unit)), enc(units.from_file(float(x), unit)),
                           units.display(float(x), unit)]
                          for x in finite if abs(x) < 1e300 for unit in units.UNITS + ('°C',)]
    out['rounds'] = [[float(x), n, round(float(x), n)] for x in finite for n in (0, 2, 8)
                     if abs(x) < 1e300] + [[2.675, 2, round(2.675, 2)], [0.125, 2, round(0.125, 2)],
                                           [0.375, 2, round(0.375, 2)], [389.60000000000002, 8, round(389.6, 8)]]
    out['plus_one'] = [[float(x), _plus_one(x)] for x in finite if abs(x) < 1e300] + \
        [[x, _plus_one(x)] for x in (1 / 3, 2 / 3, 1e-30, 4.000000000000001, 1e-17, 0.30000000000000004)]
    out['refined_temperatures'] = []
    for grid, refine in [([300, 450], {'start_C': 340, 'stop_C': 370, 'step_C': 10}),
                         ([300.0, 450.0], {'start_C': 300, 'stop_C': 400, 'step_C': 5}),
                         ([300, 400], {'start_C': 300.05, 'stop_C': 301, 'step_C': 0.1}),
                         ([0.1, 3], {'start_C': 0.1, 'stop_C': 3, 'step_C': 0.3}),
                         ([100, 580], {'start_C': 310, 'stop_C': 390, 'step_C': 7.5}),
                         ([300, 310], {'start_C': 305, 'stop_C': 305, 'step_C': 1}),
                         ([300, 310], {'start_C': 320, 'stop_C': 330, 'step_C': 1}),
                         ([300, 310], None)]:
        out['refined_temperatures'].append([grid, refine, enc(refined_temperatures(grid, refine))])
    return out


def main():
    files = sorted((ROOT / 'cases').rglob('*.toml')) + sorted((ROOT / 'tests' / 'data' / 'legacy_cases').glob('*.toml')) \
        + [ROOT / TEMPLATE]
    out = {'files': [file_dump(path) for path in files]}
    out['variants'] = []
    for label, base, replacements in VARIANTS:
        text = variant_text(base, replacements)
        out['variants'].append({'label': label, 'base': base, 'replace': [list(r) for r in replacements],
                                **text_dump(text, 'variant')})
    out['feeds'], out['legacy_feeds'] = feed_calls()
    out.update(helpers())
    # species data the tests use in place of the registry
    names = set(FEED_SPECIES)
    for path in files:
        for names_ in tomllib.loads(path.read_text(encoding='utf-8')).get('phases', {}).values():
            names.update(names_)
    known = {**MODEL.registry.all_molecular, **SOLIDS, **MODEL.registry.ion}
    out['species'] = {name: {'elements': enc(ELEMENTS_OF(name)), 'solid': name in SOLIDS}
                      for name in sorted(names) if name in known}
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(out, ensure_ascii=False, separators=(',', ':')) + '\n', encoding='utf-8')
    print(f'{OUT}: {len(files)} files, {len(VARIANTS)} variants, {OUT.stat().st_size / 1000:.0f} kB')


if __name__ == '__main__':
    main()
