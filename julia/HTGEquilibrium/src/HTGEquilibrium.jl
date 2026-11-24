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

export FormulaError, parse_formula, split_charge, molar_mass, formula_summary, hill_formula,
       parse_element_list, is_plain_formula
export WaterModel, SUPCRT92Water, IAPWS95Water, WaterProperties, water_properties, water_dict
export HKFRecord, hkf_ion_g0, hkf_apparent_G
export SUPCRT92Band, supcrt92_failed, supcrt92_failure_bands, supcrt92_band_check

include("constants.jl")
include("formulas.jl")
include("water.jl")
include("hkf.jl")
include("supcrt92_bands.jl")

end # module
