@reexport module GR4J

using ..DeepFlex
using ..DeepFlex.NamedTupleTools
"""
SoilWaterReservoir in GR4J
"""
function Surface(; name::Symbol)
    funcs = [
        DeepFlex.RainfallFlux([:prcp, :pet]),
        DeepFlex.SimpleFlux([:prcp, :pet], :pet,
            param_names=Symbol[],
            func=(i, p, sf) -> @.(sf(i[:pet] - i[:prcp]) * (i[:pet] - i[:prcp]))),
        DeepFlex.InfiltrationFlux([:rainfall])
    ]

    DeepFlex.HydroElement(
        Symbol(name, :_surf_),
        funcs=funcs
    )
end

"""
SoilWaterReservoir in GR4J
"""
function Soil(; name::Symbol)

    funcs = [
        DeepFlex.SaturationFlux([:soilwater, :infiltration], param_names=[:x1]),
        DeepFlex.EvapFlux([:soilwater, :pet], param_names=[:x1]),
        DeepFlex.PercolationFlux([:soilwater], param_names=[:x1]),
        DeepFlex.SimpleFlux([:infiltration, :percolation, :saturation], :tempflow,
            param_names=Symbol[], func=(i, p, sf) -> @.(i[:infiltration] - i[:saturation] + i[:percolation])),
        DeepFlex.SimpleFlux([:tempflow], [:slowflow, :fastflow],
            param_names=Symbol[], func=(i, p, sf) -> [i[:tempflow] .* 0.9, i[:tempflow] .* 0.1]),
    ]

    dfuncs = [
        DeepFlex.DifferFlux([:saturation], [:evap, :percolation], :soilwater)
    ]

    DeepFlex.HydroElement(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs
    )
end

function Zone(; name::Symbol)

    funcs = [
        DeepFlex.RechargeFlux([:routingstore], param_names=[:x2, :x3, :ω]),
        DeepFlex.SimpleFlux([:routingstore], :routedflow,
            param_names=[:x3, :γ],
            func=(i, p, sf) -> @.((abs(p[:x3])^(1 - p[:γ])) / (p[:γ] - 1) * (abs(i[:routingstore])^p[:γ]))),
        DeepFlex.SimpleFlux([:routedflow, :recharge, :fastflow], :flow,
            param_names=Symbol[],
            func=(i, p, sf) -> @.(i[:routedflow] + i[:recharge] + i[:fastflow]))
    ]

    dfuncs = [
        DeepFlex.DifferFlux([:slowflow, :recharge], [:routedflow], :routingstore),
    ]

    DeepFlex.HydroElement(
        Symbol(name, :_zone_),
        funcs=funcs,
        dfuncs=dfuncs
    )
end


function Route(; name::Symbol)

    funcs = [
        DeepFlex.LagFlux(:slowflow, :slowflow, lag_func=DeepFlex.uh_1_half, param_names=:x4),
        DeepFlex.LagFlux(:fastflow, :fastflow, lag_func=DeepFlex.uh_2_full, param_names=:x4),
    ]

    DeepFlex.HydroElement(
        Symbol(name, :_route_),
        funcs=funcs
    )
end


function Node(; name::Symbol)

    unit = [
        Surface(name=name),
        Soil(name=name),
        Route(name=name),
        Zone(name=name),
    ]

    DeepFlex.HydroNode(
        name,
        units=unit,
    )
end
end