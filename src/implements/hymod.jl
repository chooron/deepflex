@reexport module HyMOD

using LumpedHydro

"""
SnowWaterReservoir in HyMOD
"""
function Surface(; name::Symbol, mtk::Bool=true)
    funcs = [
        RainfallFlux([:prcp, :pet]),
        SimpleFlux([:prcp, :pet], :pet,
            param_names=Symbol[],
            func=(i, p, sf) -> @.(sf(i[:pet] - i[:prcp]) * (i[:pet] - i[:prcp]))),
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
        SimpleFlux([:saturation], :fastflow, param_names=[:a], func=(i, p, sf) -> @.(i[:saturation] * (1 - p[:a]))),
        SimpleFlux([:saturation], :slowflow, param_names=[:a], func=(i, p, sf) -> @.(i[:saturation] * p[:a]))
    ]

    dfuncs = [
        DifferFlux([:infiltration], [:evap, :saturation], :soilwater)
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
        SimpleFlux([:fr1], :qf1, param_names=[:kf], func=(i, p, sf) -> p[:kf] .* i[:fr1]),
        SimpleFlux([:fr2], :qf2, param_names=[:kf], func=(i, p, sf) -> p[:kf] .* i[:fr2]),
        SimpleFlux([:fr3], :qf3, param_names=[:kf], func=(i, p, sf) -> p[:kf] .* i[:fr3]),
        SimpleFlux([:sr], :qs, param_names=[:ks], func=(i, p, sf) -> p[:ks] .* i[:sr]),
        SimpleFlux([:qs, :qf3], :flow, param_names=Symbol[], func=(i, p, sf) -> i[:qs] .+ i[:qf3]),
    ]

    dfuncs = [
        DifferFlux([:fastflow], [:qf1], :fr1),
        DifferFlux([:qf1], [:qf2], :fr2),
        DifferFlux([:qf2], [:qf3], :fr3),
        DifferFlux([:slowflow], [:qs], :sr),
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
        surface=Surface(name=name, mtk=mtk),
        soil=Soil(name=name, mtk=mtk),
        freewater=FreeWater(name=name),
        step=step,
    )
end

function Route(; name::Symbol, mtk::Bool=true)

    funcs = [
        SimpleFlux([:flow], :flow, param_names=Symbol[], func=(i, p, sf) -> i[:flow])
    ]

    HydroElement(
        Symbol(name, :_route_),
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