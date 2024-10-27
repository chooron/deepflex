include("../src/HydroModels.jl")

using ModelingToolkit
using ComponentArrays
using Integrals
using DataInterpolations
using Plots

@variables q
@parameters lag

UnitHydroFlux = HydroModels.UnitHydroFlux
uhflux_1 = UnitHydroFlux(q, lag, :UH_1_HALF, solvetype=:DISCRETE)
uhflux_2 = UnitHydroFlux(q, lag, :UH_1_HALF, solvetype=:SPARSE)
uhflux_3 = UnitHydroFlux(q, lag, :UH_1_HALF, solvetype=:INTEGRAL)

input_mat = [1 3 4 10 14 18 23 34 23 21 17 12 6 3]
output_1 = uhflux_1(input_mat, ComponentVector(params=(lag=2.4,)))
output_2 = uhflux_2(input_mat, ComponentVector(params=(lag=2.4,)))
output_3 = uhflux_3(input_mat, ComponentVector(params=(lag=2.4,)))
uh_type_f(t, lag) =begin
    if t - lag > 0
        typeof(lag)(1)
    else
        (t / lag)^2.5
    end
end
prob = IntegralProblem(uh_type_f, (0, 10), 1.0)
sol = solve(prob, QuadGKJL())
sol.u

plot(input_mat[1, :], label="inflow")
plot!(output_1[1, :], label="uhflux_1")
plot!(output_2[1, :], label="uhflux_2")
plot!(output_3[1, :], label="uhflux_3")

