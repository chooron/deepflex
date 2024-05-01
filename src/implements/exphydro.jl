@reexport module ExpHydro

using ..DeepFlex
using ..DeepFlex.NamedTupleTools

"""
SoilWaterReservoir in Exp-Hydro
"""
function Surface(; name::Symbol, mtk::Bool=true)
    funcs = [
        DeepFlex.PetFlux([:temp, :lday]),
        DeepFlex.SnowfallFlux([:prcp, :temp], param_names=[:Tmin]),
        DeepFlex.MeltFlux([:snowwater, :temp], param_names=[:Tmax, :Df]),
        DeepFlex.RainfallFlux([:prcp, :temp], param_names=[:Tmin]),
        DeepFlex.InfiltrationFlux([:rainfall, :melt])
    ]

    dfuncs = [
        DeepFlex.DifferFlux([:snowfall], [:melt], :snowwater),
    ]

    DeepFlex.HydroElement(
        Symbol(name, :_surf_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk,
    )
end

"""
SoilWaterReservoir in Exp-Hydro
"""
function Soil(; name::Symbol, mtk::Bool=true)
    funcs = [
        DeepFlex.EvapFlux([:soilwater, :pet], param_names=[:Smax]),
        DeepFlex.BaseflowFlux([:soilwater], param_names=[:Smax, :Qmax, :f]),
        DeepFlex.SurfaceflowFlux([:soilwater], param_names=[:Smax]),
    ]

    dfuncs = [
        DeepFlex.DifferFlux([:infiltration], [:evap, :baseflow, :surfaceflow], :soilwater)
    ]

    DeepFlex.HydroElement(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end

"""
Inner Route Function in Exphydro
"""
function Zone(; name::Symbol, mtk::Bool=true)

    funcs = [
        DeepFlex.FlowFlux([:baseflow, :surfaceflow])
    ]

    DeepFlex.HydroElement(
        Symbol(name, :_zone_),
        funcs=funcs,
        mtk=mtk,
    )
end

function Route(; name::Symbol, mtk::Bool=true)

    funcs = [
        DeepFlex.SimpleFlux([:flow], :flow, param_names=Symbol[], func=(i, p, sf) -> i[:flow])
    ]

    DeepFlex.HydroElement(
        name,
        funcs=funcs,
        mtk=mtk,
    )
end

function Node(; name::Symbol, mtk::Bool=true)
    unit = [
        Surface(name=name, mtk=mtk),
        Soil(name=name, mtk=mtk),
        Zone(name=name, mtk=mtk),
        Route(name=name, mtk=mtk),
    ]

    DeepFlex.HydroNode(
        name,
        units=unit,
    )
end

end