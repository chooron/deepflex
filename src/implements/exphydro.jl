@reexport module ExpHydro

using ..LumpedHydro

"""
SoilWaterReservoir in Exp-Hydro
"""
function Surface(; name::Symbol, mtk::Bool=true)
    fluxes = [
        SimpleFlux([:temp, :lday] => [:pet]),
        SimpleFlux([:prcp, :temp] => [:snowfall], [:Tmin]),
        SimpleFlux([:snowwater, :temp] => [:melt], [:Tmax, :Df]),
        SimpleFlux([:prcp, :temp] => [:rainfall], [:Tmin]),
        SimpleFlux([:rainfall, :melt] => [:infiltration])
    ]

    dfluxes = [
        StateFlux([:snowfall] => [:melt], :snowwater, funcs=fluxes),
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
        SimpleFlux([:soilwater, :pet] => [:evap], [:Smax]),
        SimpleFlux([:soilwater] => [:baseflow], [:Smax, :Qmax, :f]),
        SimpleFlux([:soilwater] => [:surfaceflow], [:Smax]),
    ]

    dfluxes = [
        StateFlux([:infiltration] => [:evap, :baseflow, :surfaceflow], :soilwater, funcs=fluxes)
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
        SimpleFlux([:baseflow, :surfaceflow] => [:flow])
    ]

    HydroElement(
        Symbol(name, :_zone),
        funcs=fluxes,
        mtk=mtk
    )
end

function Unit(; name::Symbol, mtk::Bool=true)

    elements = [
        Surface(name=name, mtk=mtk),
        Soil(name=name, mtk=mtk),
        FreeWater(name=name, mtk=mtk)
    ]

    HydroUnit(
        name,
        elements=elements,
    )
end

end