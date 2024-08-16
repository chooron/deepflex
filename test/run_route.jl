include("../src/LumpedHydro.jl")
using ComponentArrays
using ModelingToolkit: @variables, @parameters

@variables q1 q2
@parameters x1 x2
router = LumpedHydro.UnitHydroRoute([q1, q2], [x1, x2], LumpedHydro.uh_1_half)
input_mat = hcat(Float64[1, 2, 3, 4, 2, 1, 2], Float64[1, 2, 3, 4, 2, 1, 2])'
input_vec = Float64[1, 2, 3, 4, 2, 1, 2]
# re = LumpedHydro.solve_uhfunc(input_vec, [0.2, 0.3])
# router(input_mat, ComponentVector(params=(x1=1.5, x2=3.2)))
router2 = LumpedHydro.MuskingumRoute([q1, q2])
router2(input_mat, ComponentVector(k=1.2, x=0.5, dt=1.0))
# re = LumpedHydro.solve_mskfunc(input_vec, ComponentVector(k=1.2, x=0.5, dt=1.0))

# re = LumpedHydro.solve_mskfunc(input_vec, ComponentVector(k=1.2, x=0.5, dt=1.0))


# @testset "test lag flux" begin
#     @variables a a_lag
#     @parameters lt
#     lag_flux_1 = LagFlux(a => a_lag, lt, LumpedHydro.uh_1_half)
#     @test LumpedHydro.get_input_names(lag_flux_1) == (:a,)
#     @test LumpedHydro.get_param_names(lag_flux_1) == (:lt,)
#     @test LumpedHydro.get_output_names(lag_flux_1) == (:a_lag,)
#     @test lag_flux_1(Float32[2, 3, 4, 2, 3, 1], [3.5]) ≈ [[
#         0.043634488475497855, 0.334102918508042, 1.2174967306061588,
#         2.519953682639187, 3.2301609643779736, 2.7991762465729138
#     ]]
# end
