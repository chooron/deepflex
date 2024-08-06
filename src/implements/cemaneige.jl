@reexport module Cemaneige

using ..LumpedHydro
using ..LumpedHydro: ifelse_func
using ..LumpedHydro.Symbolics: @variables
using ..LumpedHydro.ModelingToolkit: @parameters
using ..LumpedHydro.ModelingToolkit: Num

"""
SoilWaterReservoir in GR4J
"""
function SurfaceStorage(; name::Symbol, mtk::Bool=true)
    @variables snowfall = 0.0
    @variables rainfall = 0.0

    @variables prcp = 0.0
    @variables mean_temp = 0.0
    @variables max_temp = 0.0
    @variables min_temp = 0.0

    @variables prcp_ = 0.0
    @variables mean_temp_ = 0.0
    @variables max_temp_ = 0.0
    @variables min_temp_ = 0.0

    @variables melt = 0.0
    @variables solid_frac = 0.0
    @variables infiltration = 0.0
    @variables d_thermal = 0.0

    @variables snowwater = 0.0
    @variables new_snowwater = 0.0
    @variables thermal = 0.0
    @variables new_thermal = 0.0

    @parameters height = 0.0
    @parameters altitude = 0.0
    @parameters zthresh = 0.0
    @parameters snwthresh = 0.0
    @parameters CTG = 0.0
    @parameters Kf = 0.0

    funcs = [
        LumpedHydro.SimpleFlux(
            [prcp] => [prcp_], [zthresh, altitude, height],
            flux_exprs=@.[ifelse(zthresh > altitude, prcp * exp((altitude - height) * 0.0004), prcp * ifelse(height <= zthresh, exp(zthresh - height) * 0.0004, 1))]
        ),
        LumpedHydro.SimpleFlux(
            [mean_temp, min_temp, max_temp] => [mean_temp_, min_temp_, max_temp_], [altitude, height],
            flux_exprs=@.[(altitude - height) * (-0.0065) + mean_temp, (altitude - height) * (-0.0065) + min_temp, (altitude - height) * (-0.0065) + max_temp]
        ),
        LumpedHydro.SimpleFlux(
            [prcp, mean_temp, max_temp, min_temp] => [solid_frac], [altitude, zthresh],
            flux_exprs=@.[ifelse_func(zthresh - altitude) * (ifelse_func(max_temp) * ifelse_func(-min_temp) * (1.0 - max_temp / (max_temp - min_temp)) + ifelse_func(-max_temp)) +
                          ifelse_func(altitude - zthresh) * (ifelse_func(-mean_temp) + ifelse_func(3 - mean_temp) * ifelse_func(mean_temp) * (1 - (mean_temp + 1) / 4.0))]
        ),
        LumpedHydro.SimpleFlux([prcp, solid_frac] => [snowfall, rainfall], flux_exprs=@.[prcp * solid_frac, prcp * (1 - solid_frac)]),
        LumpedHydro.SimpleFlux([thermal, mean_temp] => [new_thermal], [CTG], flux_exprs=@.[min(0.0, CTG * thermal + (1 - CTG) * mean_temp)]),
        LumpedHydro.SimpleFlux([snowfall, snowwater, new_thermal, mean_temp] => [melt], [Kf, CTG, snwthresh],
            flux_exprs=@.[ifelse((new_thermal == 0) & (mean_temp > 0), min(Kf * mean_temp, snowwater + snowfall), 0.0) * (0.9 * min(1.0, (snowwater + snowfall) / snwthresh) + 0.1)]
        ),
        LumpedHydro.SimpleFlux([snowfall, melt, snowwater] => [new_snowwater], flux_exprs=@.[snowwater + snowfall - melt]),
        LumpedHydro.SimpleFlux([rainfall, melt] => [infiltration], flux_exprs=@.[rainfall + melt]),
    ]

    dfuncs = [
        LumpedHydro.StateFlux(new_snowwater => snowwater),
        LumpedHydro.StateFlux(new_thermal => thermal),
    ]

    HydroElement(
        Symbol(name, :_surface),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end
end