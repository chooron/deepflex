@reexport module GR4J

using ..LumpedHydro
using ..LumpedHydro.NamedTupleTools
"""
SoilWaterReservoir in GR4J
"""
function Surface(; name::Symbol, mtk::Bool=true)
    funcs = [
        LumpedHydro.RainfallFlux([:prcp, :pet]),
        LumpedHydro.SimpleFlux([:prcp, :pet], :pet,
            param_names=Symbol[],
            func=(i, p; kw...) -> @.(get(kw, :smooth_func, step_func)(i[:pet] - i[:prcp]) * (i[:pet] - i[:prcp]))),
        LumpedHydro.InfiltrationFlux([:rainfall])
    ]

    LumpedHydro.HydroElement(
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
        LumpedHydro.SaturationFlux([:soilwater, :infiltration], param_names=[:x1]),
        LumpedHydro.EvapFlux([:soilwater, :pet], param_names=[:x1]),
        LumpedHydro.PercolationFlux([:soilwater], param_names=[:x1]),
        LumpedHydro.SimpleFlux([:infiltration, :percolation, :saturation], :tempflow,
            param_names=Symbol[], func=(i, p; kw...) -> @.(i[:infiltration] - i[:saturation] + i[:percolation])),
        LumpedHydro.SimpleFlux([:tempflow], [:slowflow, :fastflow],
            param_names=Symbol[], func=(i, p; kw...) -> [i[:tempflow] .* 0.9, i[:tempflow] .* 0.1]),
    ]

    dfuncs = [
        LumpedHydro.DifferFlux([:saturation], [:evap, :percolation], :soilwater)
    ]

    LumpedHydro.HydroElement(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end

function Zone(; name::Symbol, mtk::Bool=true)

    funcs = [
        LumpedHydro.RechargeFlux([:routingstore], param_names=[:x2, :x3, :ω]),
        LumpedHydro.SimpleFlux([:routingstore], :routedflow,
            param_names=[:x3, :γ],
            func=(i, p; kw...) -> @.((abs(p[:x3])^(1 - p[:γ])) / (p[:γ] - 1) * (abs(i[:routingstore])^p[:γ]))),]


    dfuncs = [
        LumpedHydro.StateFlux([:slowflow, :recharge], [:routedflow], :routingstore),
    ]

    LumpedHydro.HydroElement(
        Symbol(name, :_zone_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end


function Route(; name::Symbol)

    funcs = [
        LumpedHydro.LagFlux(:slowflow, lag_func=LumpedHydro.uh_1_half, lag_time=:x4),
        LumpedHydro.LagFlux(:fastflow, lag_func=LumpedHydro.uh_2_full, lag_time=:x4),
        LumpedHydro.SimpleFlux([:routedflow, :recharge, :fastflow], :flow, param_names=Symbol[],
            func=(i, p; kw...) -> @.(i[:routedflow] + i[:recharge] + i[:fastflow]))
    ]

    LumpedHydro.HydroElement(
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

    LumpedHydro.HydroNode(
        name,
        units=namedtuple([name], [units]),
        routes=namedtuple([name], [routes]),
        step=step,
    )
end
end