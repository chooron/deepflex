@reexport module HyMOD

using ..DeepFlex
using ..DeepFlex.NamedTupleTools

"""
SnowWaterReservoir in HyMOD
"""
function Surface(; name::Symbol, mtk::Bool=true)
    funcs = [
        DeepFlex.RainfallFlux([:prcp, :pet]),
        DeepFlex.SimpleFlux([:prcp, :pet], :pet,
            param_names=Symbol[],
            func=(i, p, sf) -> @.(sf(i[:pet] - i[:prcp]) * (i[:pet] - i[:prcp]))),
        DeepFlex.InfiltrationFlux([:rainfall])
    ]

    DeepFlex.HydroElement(
        name=name,
        funcs=funcs,
        mtk=mtk,
    )
end

"""
SoilWaterReservoir in HYMOD
"""
function Soil(; name::Symbol, mtk::Bool=true)

    funcs = [
        DeepFlex.SaturationFlux([:soilwater, :infiltration], param_names=[:Smax, :b]),
        DeepFlex.EvapFlux([:soilwater, :pet], param_names=[:Smax]),
        DeepFlex.SimpleFlux([:saturation], :fastflow, param_names=[:a], func=(i, p, sf) -> @.(i[:saturation] * (1 - p[:a]))),
        DeepFlex.SimpleFlux([:saturation], :slowflow, param_names=[:a], func=(i, p, sf) -> @.(i[:saturation] * p[:a]))
    ]

    dfuncs = [
        DeepFlex.DifferFlux([:saturation], [:evap, :saturation], :soilwater)
    ]

    DeepFlex.HydroElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk,
    )
end


function Zone(; name::Symbol, mtk::Bool=true)

    funcs = [
        DeepFlex.SimpleFlux([:fr1], :qf1, param_names=[:kf], func=(i, p, sf) -> p[:kf] .* i[:fr1]),
        DeepFlex.SimpleFlux([:fr2], :qf2, param_names=[:kf], func=(i, p, sf) -> p[:kf] .* i[:fr2]),
        DeepFlex.SimpleFlux([:fr3], :qf3, param_names=[:kf], func=(i, p, sf) -> p[:kf] .* i[:fr3]),
        DeepFlex.SimpleFlux([:sr], :qs, param_names=[:ks], func=(i, p, sf) -> p[:ks] .* i[:sr]),
        DeepFlex.SimpleFlux([:qs, :qf3], :flow, param_names=Symbol[], func=(i, p, sf) -> i[:qs] .+ i[:qf3]),
    ]

    dfuncs = [
        DeepFlex.DifferFlux([:fastflow], [:qf1], :fr1),
        DeepFlex.DifferFlux([:qf1], [:qf2], :fr2),
        DeepFlex.DifferFlux([:qf2], [:qf3], :fr3),
        DeepFlex.DifferFlux([:slowflow], [:qs], :sr),
    ]

    DeepFlex.HydroElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk,
    )
end

function Route(; name::Symbol, mtk::Bool=true)

    funcs = [
        DeepFlex.SimpleFlux([:flow], :flow, param_names=Symbol[], func=(i, p, sf) -> i[:flow])
    ]

    DeepFlex.HydroElement(
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

    DeepFlex.HydroNode(
        name,
        units=namedtuple([name], [units]),
        routes=namedtuple([name], [routes]),
        step=step,
    )
end

end