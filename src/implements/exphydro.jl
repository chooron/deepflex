@reexport module ExpHydro

using ..LumpedHydro

"""
SoilWaterReservoir in Exp-Hydro
"""
function Surface(; name::Symbol, mtk::Bool=true)
    fluxes = [
        SimpleFlux([:temp, :lday], [:pet]),
        SimpleFlux([:prcp, :temp], [:snowfall], param_names=[:Tmin]),
        SimpleFlux([:snowwater, :temp], [:melt], param_names=[:Tmax, :Df]),
        SimpleFlux([:prcp, :temp], [:rainfall], param_names=[:Tmin]),
        SimpleFlux([:rainfall, :melt], [:infiltration])
    ]

    dfluxes = [
        StateFlux([:snowfall], [:melt], :snowwater, fluxes=fluxes),
    ]

    HydroElement(
        Symbol(name, :_surface),
        funcs=fluxes,
        dfuncs=dfluxes,
        mtk=mtk,
    )
end

"""
SoilWaterReservoir in Exp-Hydro
"""
function Soil(; name::Symbol, mtk::Bool=true)
    fluxes = [
        SimpleFlux([:soilwater, :pet], [:evap], param_names=[:Smax]),
        SimpleFlux([:soilwater], [:baseflow], param_names=[:Smax, :Qmax, :f]),
        SimpleFlux([:soilwater], [:surfaceflow], param_names=[:Smax]),
    ]

    dfluxes = [
        StateFlux([:infiltration], [:evap, :baseflow, :surfaceflow], :soilwater, fluxes=fluxes)
    ]

    HydroElement(
        Symbol(name, :_soil),
        funcs=fluxes,
        dfuncs=dfluxes,
        mtk=mtk
    )
end

"""
Inner Route Function in Exphydro
"""
function FreeWater(; name::Symbol, mtk=true)

    fluxes = [
        SimpleFlux([:baseflow, :surfaceflow], :totalflow)
    ]

    HydroElement(
        Symbol(name, :_zone),
        funcs=fluxes,
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
        LagFlux(:totalflow, :flow, param_names=Symbol[], func=(i, p; kw...) -> i[:totalflow])
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