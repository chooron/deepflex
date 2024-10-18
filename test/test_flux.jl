include("../src/HydroModels.jl")

using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ComponentArrays
using Lux
using StableRNGs

chain = Lux.Chain(
    Lux.Dense(2, 2),
    Lux.Dense(2, 1),
    name=:nn
)

@variables a b c d
@parameters p1 p2
tv_flux_1 = HydroModels.SimpleFlux([a, b] => [c], [p1, p2], exprs=[a * p1 + b * p2])
state_flux = HydroModels.StateFlux([a, b] => [c], d)
nn_flux =HydroModels.NeuralFlux([a, b] => [c], chain)
mr_flux = HydroModels.MuskingumRouteFlux(a)
dc_flux = HydroModels.DischargeRouteFlux(a)
uh_flux = HydroModels.UnitHydroFlux(a, p1, HydroModels.uh_1_half)