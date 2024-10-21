@reexport module ExpHydro

using ..HydroModels

"""
SoilWaterReservoir in Exp-Hydro
"""
function SurfaceStorage(; name::Symbol)
    fluxes = [
        SimpleFlux([:temp, :lday] => [:pet]),
        SimpleFlux([:prcp, :temp] => [:snowfall], [:Tmin]),
        SimpleFlux([:snowwater, :temp] => [:melt], [:Tmax, :Df]),
        SimpleFlux([:prcp, :temp] => [:rainfall], [:Tmin]),
        SimpleFlux([:rainfall, :melt] => [:infiltration])
    ]

    dfluxes = [
        StateFlux([:snowfall] => [:melt], :snowwater),
    ]

    HydroBucket(
        name=Symbol(name, :_surface),
        funcs=fluxes,
        dfuncs=dfluxes,
    )
end

"""
SoilWaterReservoir in Exp-Hydro
"""
function SoilStorage(; name::Symbol)
    fluxes = [
        SimpleFlux([:soilwater, :pet] => [:evap], [:Smax]),
        SimpleFlux([:soilwater] => [:baseflow], [:Smax, :Qmax, :f]),
        SimpleFlux([:soilwater] => [:surfaceflow], [:Smax]),
    ]

    dfluxes = [
        StateFlux([:infiltration] => [:evap, :baseflow, :surfaceflow], :soilwater)
    ]

    HydroBucket(
        name=Symbol(name, :_soil),
        funcs=fluxes,
        dfuncs=dfluxes,
    )
end

"""
Inner Route Function in Exphydro
"""
function FreeWater(; name::Symbol)

    fluxes = [
        SimpleFlux([:baseflow, :surfaceflow] => [:flow], flux_funcs=[(i, p) -> i[1] + i[2]])
    ]

    HydroBucket(
        name=Symbol(name, :_zone),
        funcs=fluxes,
    )
end

function Model(; name::Symbol)

    elements = [
        SurfaceStorage(name=name),
        SoilStorage(name=name),
        FreeWater(name=name)
    ]

    HydroModel(
        name,
        components=elements,
    )
end

end