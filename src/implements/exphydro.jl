@reexport module ExpHydro

using ..HydroModels
using ..HydroModels: @variables, @parameters
using ..HydroModels: step_func

"""
SoilWaterReservoir in Exp-Hydro
"""
function SurfaceStorage(; name::Symbol)
    @variables temp lday pet prcp snowfall rainfall snowpack melt
    @parameters Tmin Tmax Df

    fluxes = [
        SimpleFlux([temp, lday] => [pet],
            exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
        SimpleFlux([prcp, temp] => [snowfall, rainfall], [Tmin],
            exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
        SimpleFlux([snowpack, temp] => [melt], [Tmax, Df],
            exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
    ]

    dfluxes = [
        StateFlux([snowfall] => [melt], snowpack),
    ]

    HydroBucket(
        Symbol(name, :_surface),
        funcs=fluxes,
        dfuncs=dfluxes,
    )
end

"""
SoilWaterReservoir in Exp-Hydro
"""
function SoilStorage(; name::Symbol)
    @variables soilwater pet evap baseflow surfaceflow flow rainfall melt
    @parameters Smax Qmax f
    fluxes = [
        SimpleFlux([soilwater, pet] => [evap], [Smax],
            exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
        SimpleFlux([soilwater] => [baseflow], [Smax, Qmax, f],
            exprs=[step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))]),
        SimpleFlux([soilwater] => [surfaceflow], [Smax],
            exprs=[max(0.0, soilwater - Smax)]),
        SimpleFlux([baseflow, surfaceflow] => [flow],
            exprs=[baseflow + surfaceflow]),
    ]

    dfluxes = [
        StateFlux([rainfall, melt] => [evap, flow], soilwater)
    ]

    HydroBucket(
        Symbol(name, :_soil),
        funcs=fluxes,
        dfuncs=dfluxes,
    )
end

function Model(; name::Symbol)

    elements = [
        SurfaceStorage(name=name),
        SoilStorage(name=name),
    ]

    HydroModel(
        name,
        components=elements,
    )
end
end
