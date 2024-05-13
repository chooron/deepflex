@reexport module HyMOD

using ..LumpedHydro
using ..LumpedHydro.NamedTupleTools

"""
SnowWaterReservoir in HyMOD
"""
function Surface(; name::Symbol, mtk::Bool=true)
    funcs = [
        LumpedHydro.RainfallFlux([:prcp, :pet]),
        LumpedHydro.SimpleFlux([:prcp, :pet], :pet,
            param_names=Symbol[],
            func=(i, p, sf) -> @.(sf(i[:pet] - i[:prcp]) * (i[:pet] - i[:prcp]))),
        LumpedHydro.InfiltrationFlux([:rainfall])
    ]

    LumpedHydro.HydroElement(
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
        LumpedHydro.SaturationFlux([:soilwater, :infiltration], param_names=[:Smax, :b]),
        LumpedHydro.EvapFlux([:soilwater, :pet], param_names=[:Smax]),
        LumpedHydro.SimpleFlux([:saturation], :fastflow, param_names=[:a], func=(i, p, sf) -> @.(i[:saturation] * (1 - p[:a]))),
        LumpedHydro.SimpleFlux([:saturation], :slowflow, param_names=[:a], func=(i, p, sf) -> @.(i[:saturation] * p[:a]))
    ]

    dfuncs = [
        LumpedHydro.DifferFlux([:infiltration], [:evap, :saturation], :soilwater)
    ]

    LumpedHydro.HydroElement(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk,
    )
end


function Zone(; name::Symbol, mtk::Bool=true)

    funcs = [
        LumpedHydro.SimpleFlux([:fr1], :qf1, param_names=[:kf], func=(i, p, sf) -> p[:kf] .* i[:fr1]),
        LumpedHydro.SimpleFlux([:fr2], :qf2, param_names=[:kf], func=(i, p, sf) -> p[:kf] .* i[:fr2]),
        LumpedHydro.SimpleFlux([:fr3], :qf3, param_names=[:kf], func=(i, p, sf) -> p[:kf] .* i[:fr3]),
        LumpedHydro.SimpleFlux([:sr], :qs, param_names=[:ks], func=(i, p, sf) -> p[:ks] .* i[:sr]),
        LumpedHydro.SimpleFlux([:qs, :qf3], :flow, param_names=Symbol[], func=(i, p, sf) -> i[:qs] .+ i[:qf3]),
    ]

    dfuncs = [
        LumpedHydro.DifferFlux([:fastflow], [:qf1], :fr1),
        LumpedHydro.DifferFlux([:qf1], [:qf2], :fr2),
        LumpedHydro.DifferFlux([:qf2], [:qf3], :fr3),
        LumpedHydro.DifferFlux([:slowflow], [:qs], :sr),
    ]

    LumpedHydro.HydroElement(
        Symbol(name, :_zone_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk,
    )
end

function Route(; name::Symbol, mtk::Bool=true)

    funcs = [
        LumpedHydro.SimpleFlux([:flow], :flow, param_names=Symbol[], func=(i, p, sf) -> i[:flow])
    ]

    LumpedHydro.HydroElement(
        Symbol(name, :_route_),
        funcs=funcs,
        mtk=mtk
    )
end

function Node(; name::Symbol, mtk::Bool=true, step::Bool=true)
    units = [
        Surface(name=name,mtk=mtk),
        Soil(name=name,mtk=mtk),
        Zone(name=name,mtk=mtk),
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