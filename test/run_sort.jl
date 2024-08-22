# @testset "test flux sort function" begin
#     @variables a b c d e f
#     @parameters p1 p2 p3 p4

#     flux_1 = HydroModels.SimpleFlux([a, b] => [c, d], [p1, p2], exprs=[a * p1 + p2, b * p2 + p1])
#     flux_2 = HydroModels.SimpleFlux([a, c] => [e], [p3], exprs=[a * p3 + c])
#     flux_3 = HydroModels.SimpleFlux([e, d] => [f], exprs=[e + d])
#     sorted_fluxes = HydroModels.sort_fluxes([flux_2, flux_3, flux_1])
#     @test sorted_fluxes == [flux_1, flux_2, flux_3]
# end

# @testset "test buckets sort function" begin
include("../src/HydroModels.jl")
HydroBucket = HydroModels.HydroBucket
SimpleFlux = HydroModels.SimpleFlux
StateFlux = HydroModels.StateFlux
step_func = HydroModels.step_func
@variables temp lday pet prcp snowfall rainfall snowpack melt
@parameters Tmin Tmax Df
name =:test
snow_fluxes = [
    SimpleFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    SimpleFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    SimpleFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
]
snow_dfluxes = [StateFlux([snowfall] => [melt], snowpack),]
snow_bucket = HydroBucket(Symbol(name, :_surface), funcs=snow_fluxes, dfuncs=snow_dfluxes)

@variables soilwater pet evap baseflow surfaceflow flow rainfall melt
@parameters Smax Qmax f
soil_fluxes = [
    SimpleFlux([soilwater, pet] => [evap], [Smax], exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
    SimpleFlux([soilwater] => [baseflow], [Smax, Qmax, f], exprs=[step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))]),
    SimpleFlux([soilwater] => [surfaceflow], [Smax], exprs=[max(0.0, soilwater - Smax)]),
]
soil_dfluxes = [StateFlux([rainfall, melt] => [evap, flow], soilwater)]
soil_bucket = HydroBucket(Symbol(name, :_soil), funcs=soil_fluxes, dfuncs=soil_dfluxes)
flow_flux = SimpleFlux([baseflow, surfaceflow] => [flow], exprs=[baseflow + surfaceflow])
sorted_fluxes = HydroModels.sort_components([flow_flux, soil_bucket, snow_bucket])