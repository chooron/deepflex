using HydroModels
using ComponentArrays

@variables temp lday

pet_flux = HydroFlux([temp,lday]=>[pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)])

