include("../src/HydroModels.jl")
using ModelingToolkit
using ComponentArrays
using Graphs

@variables q1

k = 3.0
x = 0.2
dt = 1.0
pas = ComponentVector(params=(k=k, x=x, dt=dt,))

msk_flux = HydroModels.MuskingumRouteFlux(q1)

input = Float64[1 2 3 2 3 2 5 7 8 3 2 1]

re = msk_flux(input, pas)

# function msk_func(input_vec)
#     c0 = ((dt / k) - (2 * x)) / ((2 * (1 - x)) + (dt / k))
#     c1 = ((dt / k) + (2 * x)) / ((2 * (1 - x)) + (dt / k))
#     c2 = ((2 * (1 - x)) - (dt / k)) / ((2 * (1 - x)) + (dt / k))

#     q_out = input_vec[1]
#     q_out_vec = [q_out]
#     for i in 2:length(input_vec)
#         q_out = c0 * input_vec[i] + c1 * input_vec[i-1] + c2 * q_out
#         push!(q_out_vec, q_out)
#     end
#     q_out_vec
# end

# re2 = msk_func(input)'