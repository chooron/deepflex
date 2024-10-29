include("../src/HydroModels.jl")
using ModelingToolkit
using ComponentArrays
using Plots

@variables q

flux_1 = HydroModels.MuskingumRouteFlux(q)
flux_2 = HydroModels.RiverRouteFlux(q)

inflow = Float64[4 6 8 14 23 36 45 32 21 12 8]
timeidx = collect(1:length(inflow))

params = ComponentVector(k=1.2, x=0.1)
pas = ComponentVector(params=params)

q_out_1 = flux_1(inflow, pas, timeidx=timeidx, delta_t=1.0)
q_out_2 = flux_2(inflow, pas, timeidx=timeidx, delta_t=1.0)

plot(inflow[1, :])
plot!(q_out_1)
plot!(q_out_2)
