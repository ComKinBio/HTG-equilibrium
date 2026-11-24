"""
HTGEquilibrium: the science core of the HTG equilibrium app in Julia (Phase 2,
milestone M2). Python module names are kept (htg/constants.py -> constants.jl, ...)
so the two versions are easy to compare.
"""
module HTGEquilibrium

using Clapeyron
using ForwardDiff
using Libdl
using OrderedCollections
using Printf
using TOML

export FormulaError, parse_formula, split_charge, molar_mass, formula_summary, hill_formula,
       parse_element_list, is_plain_formula
export WaterModel, SUPCRT92Water, IAPWS95Water, WaterProperties, water_properties, water_dict
export HKFRecord, hkf_ion_g0, hkf_apparent_G
export SUPCRT92Band, supcrt92_failed, supcrt92_failure_bands, supcrt92_band_check
export ValueError, pyrepr, parse_toml, read_toml, Units
export Salt, Feed, Case, CaseError, LegacyConversionError, load_case, case_from_dict, validate,
       temperatures, pressures, grid, chains, at, is_single_point, resolved, total_water_g_per_g_dry,
       dry_element_wt_percent, carbon_bound_messages, effective_carbon_bounds, check_carbon_bounds,
       legacy_carbon_bounds, legacy_to_v2, is_v2, fluid_phase_count, starts_from_multistart,
       BOUND_KINDS, PHASE_LABELS
export FeedError, feed_v2_to_model, feed_summary, feed_fractions, molar_mass_function,
       species_elements_function, elemental_masses, ash_percent, scaled_elements, legacy_feed_to_v2

include("constants.jl")
include("formulas.jl")
include("water.jl")
include("hkf.jl")
include("supcrt92_bands.jl")
include("pycompat.jl")
include("toml_ordered.jl")
include("units.jl")
include("case.jl")
include("feed_v2.jl")

end # module
