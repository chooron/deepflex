@reexport module ExpHydro

using LumpedHydro

"""
SoilWaterReservoir in Exp-Hydro
"""
function Surface(; name::Symbol, mtk::Bool=true)
    funcs = [
        PetFlux([:temp, :lday]),
        SnowfallFlux([:prcp, :temp], param_names=[:Tmin]),
        MeltFlux([:snowwater, :temp], param_names=[:Tmax, :Df]),
        RainfallFlux([:prcp, :temp], param_names=[:Tmin]),
        InfiltrationFlux([:rainfall, :melt])
    ]

    dfuncs = [
        StateFlux([:snowfall], [:melt], :snowwater),
    ]

    HydroElement(
        Symbol(name, :_surface_),
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
        EvapFlux([:soilwater, :pet], param_names=[:Smax]),
        BaseflowFlux([:soilwater], param_names=[:Smax, :Qmax, :f]),
        SurfaceflowFlux([:soilwater], param_names=[:Smax]),
    ]

    dfuncs = [
        StateFlux([:infiltration], [:evap, :baseflow, :surfaceflow], :soilwater)
    ]

    HydroElement(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end

"""
Inner Route Function in Exphydro
"""
function FreeWater(; name::Symbol)

    funcs = [
        FlowFlux([:baseflow, :surfaceflow], :totalflow)
    ]

    HydroElement(
        Symbol(name, :_zone_),
        funcs=funcs,
        mtk=false
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

function Route(; name::Symbol)
    
    funcs = [
        SimpleFlux(:totalflow, :flow, param_names=Symbol[], func=(i, p; kw...) -> i[:totalflow])
    ]

    HydroElement(
        name,
        funcs=funcs,
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