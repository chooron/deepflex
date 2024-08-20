include("../src/HydroModels.jl")
using ComponentArrays
using ModelingToolkit: @variables, @parameters
using Test


@variables q1 q2
@parameters x1 x2
router = HydroModels.UnitHydroRouteFlux(q1, x1, HydroModels.uh_1_half)
input = [2 3 4 2 3 1]
router(input, ComponentVector(x1=3.5))
# input_mat = hcat(Float64[1, 2, 3, 4, 2, 1, 2], Float64[1, 2, 3, 4, 2, 1, 2])'
# input_vec = Float64[1, 2, 3, 4, 2, 1, 2]


# re = LumpedHydro.solve_uhfunc(input_vec, [0.2, 0.3])
# router(input_mat, ComponentVector(params=(x1=1.5, x2=3.2)))
# router2 = LumpedHydro.MuskingumRoute([q1, q2])
# router2(input_mat, ComponentVector(k=1.2, x=0.5, dt=1.0))
# re = LumpedHydro.solve_mskfunc(input_vec, ComponentVector(k=1.2, x=0.5, dt=1.0))

# re = LumpedHydro.solve_mskfunc(input_vec, ComponentVector(k=1.2, x=0.5, dt=1.0))


@testset "test lag flux" begin
    @variables q1
    @parameters x1
    router = HydroModels.UnitHydroRouteFlux(q1, x1, HydroModels.uh_1_half)
    @test HydroModels.get_input_names(router) == [:q1]
    @test HydroModels.get_param_names(router) == [:x1]
    @test HydroModels.get_output_names(router) == [:q1_routed]
    @test router(Float32[2 3 4 2 3 1], ComponentVector(x1=3.5)) ≈ [0.043634488475497855 0.334102918508042 1.2174967306061588 2.519953682639187 3.2301609643779736 2.7991762465729138]
end
