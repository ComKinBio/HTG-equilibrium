# Standard-state Gibbs energies g°(T, P) [J/mol] (Python: htg/standard_states.py).
#
# All sources share one reference: Hf(298.15 K) - T S, as used by NASA-7.
#
#     molecules       NASA-7 polynomials (pressure independent, ideal-gas standard state)
#     NASA solids     NASA-7 polynomials
#     SUPCRT minerals Maier-Kelley Cp = a + bT + c/T^2 plus a constant-volume term Vr (P - Pr)
#     aqueous ions    revised HKF (hkf.jl), shifted from the apparent Gf(298.15 K) anchor
#                     to the NASA reference
#
# Water properties (Python `supcrt_water`) are `water_properties(model, T, P)`
# in water.jl. Ions take the water model as their first argument: one model
# (SUPCRT92Water() by default, or IAPWS95Water()) is used throughout, never mixed.

"""
    nasa_g0(registry, nasa_db, species, T) -> Vector{Float64}

g° [J/mol] of fluid molecules (or NASA solids) from their NASA-7 records at T [K].
"""
function nasa_g0(registry::Registry, nasa_db, species, T)
    mu0 = zeros(Float64, length(species))
    for (i, name) in enumerate(species)
        record = nasa_name(registry, name)
        _has(nasa_db, record) ||
            throw(KeyError("NASA species '$(py_str(record))' for '$name' is not present in nasa_db."))
        mu0[i] = nasa7_properties(nasa_db, record, T)[4]
    end
    return mu0
end

"""
    supcrt_mineral_properties(params, name, T, P) -> (G, H, S, Cp, V_m3_per_mol)

Mineral G, H [J/mol], S, Cp [J/mol/K] and V [m3/mol] at T [K], P [Pa] from the
single-range Maier-Kelley record `params` (Gf, Hf, Sr, Vr, a, b, c, Tmax in
supcrt92.xml units):

    H(T) = Hf(Tr) + int Cp dT,  S(T) = Sr + int Cp/T dT,  G = H - T S + Vr (P - Pr).

Tmax is the equation limit, not a hydrate stability temperature.
"""
function supcrt_mineral_properties(params, name, T, P)
    t, tr, p = Float64(T), T_REF, Float64(P)
    if !isfinite(t) || !(tr <= t <= params["Tmax"])
        throw(ArgumentError("$name: T=$(_g(t)) K outside supported [$(_g(tr)), $(_g(params["Tmax"]))] K"))
    end
    (!isfinite(p) || p <= 0) && throw(ArgumentError("Invalid mineral pressure"))
    a, b, c = (params[k] * CAL_TO_J for k in ("a", "b", "c"))
    dh = a * (t - tr) + b * (t * t - tr * tr) / 2 - c * (1 / t - 1 / tr)
    ds = a * log(t / tr) + b * (t - tr) - c * (pypow(t, -2) - pypow(tr, -2)) / 2
    h, entropy = params["Hf"] * CAL_TO_J + dh, params["Sr"] * CAL_TO_J + ds
    pv = params["Vr"] * 1e-6 * (p - P_REF)
    return (G=h - t * entropy + pv, H=h + pv, S=entropy, Cp=a + b * t + c / pypow(t, 2),
            V_m3_per_mol=params["Vr"] * 1e-6)
end

"""
    solid_g0(registry, nasa_db, species, T, P) -> Vector{Float64}

g° [J/mol] of pure solids at T [K], P [Pa]: SUPCRT minerals or NASA-7 solids.
"""
function solid_g0(registry::Registry, nasa_db, species, T, P)
    values = Float64[]
    for name in species
        d = registry.solid[name]
        value = if get(d, "adapter", nothing) == "supcrt_mineral"
            supcrt_mineral_properties(d["supcrt_mineral"], name, T, P).G
        else
            nasa7_properties(nasa_db, nasa_name(registry, name), T)[4]
        end
        push!(values, value)
    end
    return values
end

"""
    hkf_records(hkf_table, ions) -> Vector{HKFRecord}

The HKF records of `ions` from the rows of `registry.ion_hkf_table`.
"""
function hkf_records(hkf_table::AbstractVector{<:AbstractDict}, ions)
    rows = Dict(String(row["Species"]) => row for row in hkf_table)
    return [haskey(rows, s) ? HKFRecord(rows[s]) : throw(KeyError(s)) for s in ions]
end

"""
    hkf_ion_g0(water_model, hkf_table, ions, T, P) -> Vector{Float64}

g° [J/mol] of aqueous ions at T [K], P [Pa] from their HKF records (rows of
`registry.ion_hkf_table`), on the NASA (Hf - T S) reference, with all water
properties from `water_model`. Python: `standard_states.hkf_ion_g0(hkf_table, ions, T, P)`.
"""
hkf_ion_g0(model::WaterModel, hkf_table::AbstractVector{<:AbstractDict}, ions, T::Real, P::Real) =
    hkf_ion_g0(model, hkf_records(hkf_table, ions), T, P)
