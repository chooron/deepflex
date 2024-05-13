@reexport module ExpHydro

using ..LumpedHydro
using ..LumpedHydro.NamedTupleTools

"""
SoilWaterReservoir in Exp-Hydro
"""
function Surface(; name::Symbol, mtk::Bool=true)
    funcs = [
        LumpedHydro.PetFlux([:temp, :lday]),
        LumpedHydro.SnowfallFlux([:prcp, :temp], param_names=[:Tmin]),
        LumpedHydro.MeltFlux([:snowwater, :temp], param_names=[:Tmax, :Df]),
        LumpedHydro.RainfallFlux([:prcp, :temp], param_names=[:Tmin]),
        LumpedHydro.InfiltrationFlux([:rainfall, :melt])
    ]

    dfuncs = [
        LumpedHydro.StateFlux([:snowfall], [:melt], :snowwater),
    ]

    LumpedHydro.HydroElement(
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
        LumpedHydro.EvapFlux([:soilwater, :pet], param_names=[:Smax]),
        LumpedHydro.BaseflowFlux([:soilwater], param_names=[:Smax, :Qmax, :f]),
        LumpedHydro.SurfaceflowFlux([:soilwater], param_names=[:Smax]),
    ]

    dfuncs = [
        LumpedHydro.StateFlux([:infiltration], [:evap, :baseflow, :surfaceflow], :soilwater)
    ]

    LumpedHydro.HydroElement(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end

"""
Inner Route Function in Exphydro
"""
function Zone(; name::Symbol)

    funcs = [
        LumpedHydro.FlowFlux([:baseflow, :surfaceflow], :totalflow)
    ]

    LumpedHydro.HydroElement(
        Symbol(name, :_zone_),
        funcs=funcs,
    )
end

function Route(; name::Symbol)

    funcs = [
        LumpedHydro.SimpleFlux(:totalflow, :flow, param_names=Symbol[], func=(i, p; kw...) -> i[:totalflow])
    ]

    LumpedHydro.HydroElement(
        name,
        funcs=funcs,
    )
end

function Node(; name::Symbol, mtk::Bool=true, step::Bool=true)
    layers = [
        Surface(name=name, mtk=mtk),
        Soil(name=name, mtk=mtk),
        Zone(name=name),
    ]

    routes = Route(name=name)

    LumpedHydro.HydroNode(
        name,
        layers=namedtuple([name], [layers]),
        routes=namedtuple([name], [routes]),
        step=step,
    )
end

end