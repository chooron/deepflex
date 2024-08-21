include("../src/HydroModels.jl")
using ComponentArrays
using ModelingToolkit: @variables, @parameters
using Test

@testset "test unit hydro flux" begin
    @variables q1
    @parameters x1
    router = HydroModels.UnitHydroRouteFlux(q1, x1, HydroModels.uh_1_half)
    @test HydroModels.get_input_names(router) == [:q1]
    @test HydroModels.get_param_names(router) == [:x1]
    @test HydroModels.get_output_names(router) == [:q1_routed]
    @test router(Float32[2 3 4 2 3 1], ComponentVector(x1=3.5)) ≈ [0.043634488475497855 0.334102918508042 1.2174967306061588 2.519953682639187 3.2301609643779736 2.7991762465729138]
end
