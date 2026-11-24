# Compute option-2 water properties and ion G on the pychnosz grid.
#   julia --project=julia/spikes/hkf_water julia/spikes/hkf_water/compare.jl
# Reads pychnosz_reference.json (dump_pychnosz.py), writes julia_results.json.
using JSON

include(joinpath(@__DIR__, "hkf_water.jl"))
using .HKFWater

const HERE = @__DIR__

ws_dict(w) = Dict(f => getfield(w, f) for f in fieldnames(HKFWater.WaterState))

function main()
    ref = JSON.parsefile(joinpath(HERE, "pychnosz_reference.json"))
    ions = String.(ref["ions"])
    recs = Dict(k => HKFWater.hkf_records_from_dict(v) for (k, v) in ref["hkf_records"])
    trpr = ref["TrPr"]
    # Reference state of water at Tr, Pr.
    wr_iapws = iapws95_water(298.15, 1.0)
    # SUPCRT92 values at Tr, Pr (only epsilon and Y enter G).
    wr_supcrt = HKFWater.WaterState(trpr["rho"], trpr["epsilon"], trpr["QBorn"], trpr["XBorn"],
                                    trpr["YBorn"], NaN, NaN, NaN)
    out = Any[]
    t0 = time()
    for p in ref["points"]
        T = p["T_C"] + 273.15
        P = Float64(p["P_bar"])
        sup = p["water"]
        # (a) implementation check: JN91 + Born from SUPCRT92's own density/alpha/beta/daldT
        wcheck = water_from_density(T, sup["rho"], sup["alpha"], sup["daldT"], sup["beta"])
        # (b) option 2: IAPWS-95 density from Clapeyron
        w2 = iapws95_water(T, P)
        row = Dict{String,Any}("T_C" => p["T_C"], "P_bar" => P,
                               "water_check" => ws_dict(wcheck), "water_iapws" => ws_dict(w2))
        mu = Dict{String,Any}()
        for ion in ions
            r = recs[ion]
            mu[ion] = Dict(
                # SUPCRT92 density/alpha/beta/daldT, Julia JN91 + g-function + HKF
                "check" => HKFWater.app_mu0(r, hkf_G(r, T, P, wcheck, wr_supcrt)),
                # option 2, everything from IAPWS-95 + JN91
                "iapws" => HKFWater.app_mu0(r, hkf_G(r, T, P, w2, wr_iapws)),
                # option 2 with SUPCRT92's Tr,Pr constants (epsilon_r, Y_r) kept
                "iapws_supref" => HKFWater.app_mu0(r, hkf_G(r, T, P, w2, wr_supcrt)))
        end
        row["mu0_ion"] = mu
        push!(out, row)
    end
    elapsed = time() - t0
    result = Dict("TrPr_iapws" => ws_dict(wr_iapws), "points" => out,
                  "seconds" => elapsed, "n" => length(out))
    open(joinpath(HERE, "julia_results.json"), "w") do io
        JSON.print(io, result, 1)
    end
    println("$(length(out)) points in $(round(elapsed, digits=2)) s; Tr,Pr IAPWS: ", wr_iapws)
end

main()
