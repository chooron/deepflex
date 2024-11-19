include("../src/HydroModels.jl")

using ModelingToolkit
using ComponentArrays
using Integrals
using DataInterpolations
using Plots

@variables q q_routed
@parameters lag

UnitHydroFlux = HydroModels.UnitHydroFlux
uhfunc = HydroModels.UHFunction(:UH_1_HALF)
uhflux_1 = UnitHydroFlux(q, q_routed, lag, uhfunc=uhfunc, solvetype=:DISCRETE)
uhflux_2 = UnitHydroFlux(q, q_routed, lag, uhfunc=uhfunc, solvetype=:SPARSE)
uhflux_3 = UnitHydroFlux(q, q_routed, lag, uhfunc=uhfunc, solvetype=:INTEGRAL)

input_mat = Float32[2 3 4 2 3 1]

output_1 = uhflux_1(input_mat, ComponentVector(params=(lag=3.5,)))
output_2 = uhflux_2(input_mat, ComponentVector(params=(lag=3.5,)))
# output_3 = uhflux_3(input_mat, ComponentVector(params=(lag=2.4,)))
# uh_type_f(t, lag) = begin
#     if t - lag > 0
#         typeof(lag)(1)
#     else
#         (t / lag)^2.5
#     end
# end
# prob = IntegralProblem(uh_type_f, (0, 10), 1.0)
# sol = solve(prob, QuadGKJL())
# sol.u

# plot(input_mat[1, :], label="inflow")
# plot!(output_1[1, :], label="uhflux_1")
# plot!(output_2[1, :], label="uhflux_2")
# plot!(output_3[1, :], label="uhflux_3")

