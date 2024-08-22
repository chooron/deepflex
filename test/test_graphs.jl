include("../src/HydroModels.jl")
using ModelingToolkit

@variables a b c d e f g
@parameters p1 p2 p3 p4

flux_1 = HydroModels.SimpleFlux([a, b] => [c, d], [p1, p2], exprs=[a * p1 + p2, b * p2 + p1])
flux_2 = HydroModels.SimpleFlux([a, c] => [e], [p3], exprs=[a * p3 + c])
flux_3 = HydroModels.SimpleFlux([e, d] => [f], exprs=[e + d])
flux_4 = HydroModels.SimpleFlux([a, b] => [g], exprs=[a + b])
sorted_fluxes = HydroModels.sort_components([flux_2, flux_4, flux_3, flux_1])