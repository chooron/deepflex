include("../src/HydroModels.jl")
using ModelingToolkit
using ComponentArrays
using Test
@variables q1
@parameters x1

router = HydroModels.UnitHydroFlux(q1, x1, HydroModels.uh_1_half, solvetype=:unithydro1)
HydroModels.get_input_names(router) == [:q1]
HydroModels.get_param_names(router) == [:x1]
HydroModels.get_output_names(router) == [:q1_routed]
router(Float32[2 2 3 2 3 4 1 2], ComponentVector(params=(x1=2.39,)))