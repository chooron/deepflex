@reexport module HyMOD

using ..LumpedHydro

"""
SnowWaterReservoir in HyMOD
"""
function Surface(; name::Symbol, mtk::Bool=true)
    funcs = [
        RainfallFlux([:prcp, :pet]),
        InfiltrationFlux([:rainfall])
    ]

    HydroElement(
        Symbol(name, :_surface_),
        funcs=funcs,
        mtk=mtk,
    )
end

"""
SoilWaterReservoir in HYMOD
"""
function Soil(; name::Symbol, mtk::Bool=true)

    funcs = [
        SaturationFlux([:soilwater, :infiltration], param_names=[:Smax, :b]),
        EvapFlux([:soilwater, :pet], param_names=[:Smax]),
        SimpleFlux([:saturation], :fastflow, param_names=[:a], func=(i, p; kw...) -> @.(i[:saturation] * (1 - p[:a]))),
        SimpleFlux([:saturation], :slowflow, param_names=[:a], func=(i, p; kw...) -> @.(i[:saturation] * p[:a]))
    ]

    dfuncs = [
        StateFlux([:infiltration], [:evap, :saturation], :soilwater)
    ]

    HydroElement(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk,
    )
end


function FreeWater(; name::Symbol, mtk::Bool=true)

    funcs = [
        SimpleFlux([:fr1], :qf1, param_names=[:kf], func=(i, p; kw...) -> p[:kf] .* i[:fr1]),
        SimpleFlux([:fr2], :qf2, param_names=[:kf], func=(i, p; kw...) -> p[:kf] .* i[:fr2]),
        SimpleFlux([:fr3], :qf3, param_names=[:kf], func=(i, p; kw...) -> p[:kf] .* i[:fr3]),
        SimpleFlux([:sr], :qs, param_names=[:ks], func=(i, p; kw...) -> p[:ks] .* i[:sr]),
        SimpleFlux([:qs, :qf3], :totalflow, param_names=Symbol[], func=(i, p; kw...) -> i[:qs] .+ i[:qf3]),
    ]

    dfuncs = [
        StateFlux([:fastflow], [:qf1], :fr1),
        StateFlux([:qf1], [:qf2], :fr2),
        StateFlux([:qf2], [:qf3], :fr3),
        StateFlux([:slowflow], [:qs], :sr),
    ]

    HydroElement(
        Symbol(name, :_zone_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk,
    )
end

function Unit(; name::Symbol, mtk::Bool=true, step::Bool=true)
    HydroUnit(
        name,
        elements=[Surface(name=name, mtk=mtk), Soil(name=name, mtk=mtk), FreeWater(name=name, mtk=mtk)],
        step=step,
    )
end

function Route(; name::Symbol, mtk::Bool=true)

    funcs = [
        SimpleFlux([:totalflow], :flow, param_names=Symbol[], func=(i, p; kw...) -> i[:totalflow])
    ]

    HydroElement(
        name,
        funcs=funcs,
        mtk=mtk
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