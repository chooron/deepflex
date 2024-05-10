@reexport module GR4J

using ..DeepFlex
using ..DeepFlex.NamedTupleTools
"""
SoilWaterReservoir in GR4J
"""
function Surface(; name::Symbol, mtk::Bool=true)
    funcs = [
        DeepFlex.RainfallFlux([:prcp, :pet]),
        DeepFlex.SimpleFlux([:prcp, :pet], :pet,
            param_names=Symbol[],
            func=(i, p; kw...) -> @.(get(kw, :smooth_func, step_func)(i[:pet] - i[:prcp]) * (i[:pet] - i[:prcp]))),
        DeepFlex.InfiltrationFlux([:rainfall])
    ]

    DeepFlex.HydroElement(
        Symbol(name, :_surf_),
        funcs=funcs,
        mtk=mtk
    )
end

"""
SoilWaterReservoir in GR4J
"""
function Soil(; name::Symbol, mtk::Bool=true)

    funcs = [
        DeepFlex.SaturationFlux([:soilwater, :infiltration], param_names=[:x1]),
        DeepFlex.EvapFlux([:soilwater, :pet], param_names=[:x1]),
        DeepFlex.PercolationFlux([:soilwater], param_names=[:x1]),
        DeepFlex.SimpleFlux([:infiltration, :percolation, :saturation], :tempflow,
            param_names=Symbol[], func=(i, p; kw...) -> @.(i[:infiltration] - i[:saturation] + i[:percolation])),
        DeepFlex.SimpleFlux([:tempflow], [:slowflow, :fastflow],
            param_names=Symbol[], func=(i, p; kw...) -> [i[:tempflow] .* 0.9, i[:tempflow] .* 0.1]),
    ]

    dfuncs = [
        DeepFlex.DifferFlux([:saturation], [:evap, :percolation], :soilwater)
    ]

    DeepFlex.HydroElement(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end

function Zone(; name::Symbol, mtk::Bool=true)

    funcs = [
        DeepFlex.RechargeFlux([:routingstore], param_names=[:x2, :x3, :ω]),
        DeepFlex.SimpleFlux([:routingstore], :routedflow,
            param_names=[:x3, :γ],
            func=(i, p; kw...) -> @.((abs(p[:x3])^(1 - p[:γ])) / (p[:γ] - 1) * (abs(i[:routingstore])^p[:γ]))),]


    dfuncs = [
        DeepFlex.StateFlux([:slowflow, :recharge], [:routedflow], :routingstore),
    ]

    DeepFlex.HydroElement(
        Symbol(name, :_zone_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end


function Route(; name::Symbol)

    funcs = [
        DeepFlex.LagFlux(:slowflow, lag_func=DeepFlex.uh_1_half, lag_time=:x4),
        DeepFlex.LagFlux(:fastflow, lag_func=DeepFlex.uh_2_full, lag_time=:x4),
        DeepFlex.SimpleFlux([:routedflow, :recharge, :fastflow], :flow, param_names=Symbol[],
            func=(i, p; kw...) -> @.(i[:routedflow] + i[:recharge] + i[:fastflow]))
    ]

    DeepFlex.HydroElement(
        Symbol(name, :_route_),
        funcs=funcs,
        mtk=false
    )
end


function Node(; name::Symbol, mtk::Bool=true, step::Bool=true)

    units = [
        Surface(name=name, mtk=mtk),
        Soil(name=name, mtk=mtk),
        Zone(name=name, mtk=mtk),
    ]

    routes = Route(name=name)

    DeepFlex.HydroNode(
        name,
        units=namedtuple([name], [units]),
        routes=namedtuple([name], [routes]),
        step=step,
    )
end
end