include("../src/HydroModels.jl")

using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ComponentArrays

@variables a b c
@parameters p1 p2
tv_flux_1 = HydroModels.TimeVaryingFlux([a, b], [c], [p1, p2], exprs=[a * p1 + b * t * p2])
tv_flux_1.func([2.0, 3.0], [1.0,2.0], [2.0])

tv_flux_1([2.0, 3.0], ComponentVector(params=(p1=1.0, p2=2.0)), timeidx=2)