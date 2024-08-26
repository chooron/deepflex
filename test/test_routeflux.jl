include("../src/HydroModels.jl")
using ComponentArrays
using ModelingToolkit
using Test

@variables q1
@parameters x1
router = HydroModels.UnitHydroRouteFlux(q1, x1, HydroModels.uh_1_half)
@test HydroModels.get_input_names(router) == [:q1]
@test HydroModels.get_param_names(router) == [:x1]
@test HydroModels.get_output_names(router) == [:q1_routed]

println(router(Float32[1 2 3 2 3 4 1 2], ComponentVector(params=(x1=2.39,))))
# [0.043634488475497855 0.334102918508042 1.2174967306061588 2.519953682639187 3.2301609643779736 2.7991762465729138]
# flow = [2, 3, 4, 5]
# weight = [2, 2, 3]
# uh_pairs = [0 => 2 .* weight, -1 => 3 .* weight, -2 => 4 .* weight, -3 => 5 .* weight]
# uh_pairs = [-(i - 1) => w .* flow for (i, w) in enumerate(weight)]
# mat = spdiagm(uh_pairs...)
# sum(mat, dims=2)