"""
HTGEquilibrium: the science core of the HTG equilibrium app in Julia (Phase 2,
milestone M2). For now it contains water properties, HKF ion standard states
and the SUPCRT92 failure-band check.
"""
module HTGEquilibrium

using Clapeyron
using ForwardDiff
using Libdl

export WaterModel, SUPCRT92Water, IAPWS95Water, WaterProperties, water_properties, water_dict
export HKFRecord, hkf_ion_g0, hkf_apparent_G
export SUPCRT92Band, supcrt92_failed, supcrt92_failure_bands, supcrt92_band_check

include("water.jl")
include("hkf.jl")
include("supcrt92_bands.jl")

end # module
