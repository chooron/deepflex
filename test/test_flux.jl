include("../src/HydroModels.jl")

using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ComponentArrays
using Lux
using StableRNGs
using Graphs

chain = Lux.Chain(
    Lux.Dense(2, 2),
    Lux.Dense(2, 1),
    name=:nn
)

@variables a b c d e f 
@parameters p1 p2 p3
tv_flux_1 = HydroModels.SimpleFlux([a, b] => [c], [p1, p2], exprs=[a * p1 + b * p2])
tv_flux_2 = HydroModels.SimpleFlux([c] => [e], [p3], exprs=[c * p3])
state_flux_1 = HydroModels.StateFlux([a, b] => [c], d)
state_flux_2 = HydroModels.StateFlux([a, b] => [e], f)
nn_flux =HydroModels.NeuralFlux([a, b] => [c], chain)
mr_flux = HydroModels.MuskingumRouteFlux(a)
dc_flux = HydroModels.DischargeRouteFlux(a)
uh_flux = HydroModels.UnitHydroFlux(a, p1, HydroModels.uh_1_half)
bucket = HydroModels.HydroBucket(name=:test, funcs=[tv_flux_1, tv_flux_2], dfuncs=[state_flux_1, state_flux_2])

# flwdir = rand(3, 3)
# positions = [(1, 1), (2, 2), (3, 3)]
# subareas = [1.0, 2.0, 3.0]
# route = HydroModels.GridRoute(:gridmr, rfunc=mr_flux, flwdir=flwdir, positions=positions, subareas=subareas)
# route = HydroModels.VectorRoute(:gridmr, rfunc=mr_flux, network=network, subareas=subareasv2)

# model = HydroModels.HydroModel(:test, components=[bucket, route])
