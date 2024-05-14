@reexport module GR4J

using ..LumpedHydro

"""
SoilWaterReservoir in GR4J
"""
function Surface(; name::Symbol, mtk::Bool=true)

    funcs = [
        RainfallFlux([:prcp, :pet]),
        InfiltrationFlux([:rainfall])
    ]

    HydroElement(
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
        SaturationFlux([:soilwater, :infiltration], param_names=[:x1]),
        EvapFlux([:soilwater, :prcp, :pet], param_names=[:x1]),
        PercolationFlux([:soilwater], param_names=[:x1]),
        SimpleFlux([:infiltration, :percolation, :saturation], :tempflow,
            param_names=Symbol[], func=(i, p; kw...) -> @.(i[:infiltration] - i[:saturation] + i[:percolation])),
        SimpleFlux([:tempflow], [:slowflow, :fastflow],
            param_names=Symbol[], func=(i, p; kw...) -> [i[:tempflow] .* 0.9, i[:tempflow] .* 0.1]),
    ]

    dfuncs = [
        StateFlux([:saturation], [:evap, :percolation], :soilwater)
    ]

    HydroElement(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end

function FreeWater(; name::Symbol, mtk::Bool=true)

    funcs = [
        RechargeFlux([:routingstore], param_names=[:x2, :x3, :ω]),
        SimpleFlux([:routingstore], :routedflow,
            param_names=[:x3, :γ],
            func=(i, p; kw...) -> @.((abs(p[:x3])^(1 - p[:γ])) / (p[:γ] - 1) * (abs(i[:routingstore])^p[:γ]))),]


    dfuncs = [
        StateFlux([:slowflow, :recharge], [:routedflow], :routingstore),
    ]

    HydroElement(
        Symbol(name, :_zone_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end

function Unit(; name::Symbol, mtk::Bool=true, step::Bool=true)
    HydroUnit(
        name,
        elements=[Surface(name=name, mtk=mtk), Soil(name=name, mtk=mtk), FreeWater(name=name, mtk=mtk)],
        step=step,
    )
end

function Route(; name::Symbol)

    funcs = [
        LagFlux(:slowflow, :slowflow_lag, lag_func=LumpedHydro.uh_1_half, lag_time=:x4),
        LagFlux(:fastflow, :fastflow_lag, lag_func=LumpedHydro.uh_2_full, lag_time=:x4),
        SimpleFlux([:routedflow, :recharge, :fastflow_lag], :flow, param_names=Symbol[],
            func=(i, p; kw...) -> @.(i[:routedflow] + i[:recharge] + i[:fastflow_lag]))
    ]

    HydroElement(
        name,
        funcs=funcs,
        mtk=false
    )
end

function Node(; name::Symbol, mtk::Bool=true, step::Bool=true)
    HydroNode(
        name,
        units=[Unit(name=name, mtk=mtk, step=step)],
        routes=[Route(name=name)],
    )
end

end