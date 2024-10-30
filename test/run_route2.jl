include("../src/HydroModels.jl")
using ModelingToolkit
using ComponentArrays
using Graphs
using Plots

solver_1 = HydroModels.DiscreteSolver()
solver_2 = HydroModels.ODESolver()


@variables q
rflux_1 = HydroModels.MuskingumRouteFlux(q)

network = DiGraph(9)
add_edge!(network, 1, 2)
add_edge!(network, 2, 5)
add_edge!(network, 3, 5)
add_edge!(network, 4, 5)
add_edge!(network, 5, 8)
add_edge!(network, 6, 9)
add_edge!(network, 7, 8)
add_edge!(network, 8, 9)

ndtypes = [:ntype1, :ntype2, :ntype3]
params = ComponentVector(NamedTuple{Tuple(ndtypes)}([(k=0.3, x=0.2) for _ in eachindex(ndtypes)]))
initstates = ComponentVector(NamedTuple{Tuple(ndtypes)}([(s_river=0.0,) for _ in eachindex(ndtypes)]))
pas = ComponentVector(; params, initstates)
vroute = HydroModels.VectorRoute(rfunc=rflux_1, network=network, subareas=1.0)
# 24 * 3600 / (10.0 * 1e6) * 1e3
input_mat = Float64[4 6 8 12 34 45 56 69 49 32 23 22 12 10 6 7 5 4 3 2 3 2 3 4 2 ]
input_arr = permutedims(repeat(input_mat, outer=[1, 1, 9]), (1, 3, 2))
timeidx = collect(1:size(input_arr)[3])
ptypes = [:ntype1, :ntype2, :ntype3, :ntype2, :ntype1, :ntype2, :ntype3, :ntype1, :ntype3]
sol_2 = vroute(input_arr, pas, config=(timeidx=(timeidx .-1) .*24, ptypes=ptypes, solver=solver_1, delta_t=24.0))

plot(sol_2[1, :, :]')

# function run_grid_route()
#     @variables q
#     rflux = HydroModels.DischargeRouteFlux(q)
#     flwdir = [1 4 8; 1 4 4; 1 1 2]
#     positions = [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]

#     params = ComponentVector(NamedTuple{Tuple(ndtypes)}([(lag=0.2,) for _ in eachindex(ndtypes)]))
#     initstates = ComponentVector(NamedTuple{Tuple(ndtypes)}([(s_river=0.1,) for _ in eachindex(ndtypes)]))
#     pas = ComponentVector(; params, initstates)
#     route = HydroModels.GridRoute(rfunc=rflux, flwdir=flwdir, positions=positions, subareas=10.0)

#     input_arr = ones(1, 9, 20)
#     timeidx = collect(1:20)
#     node_types = [:ntype1, :ntype2, :ntype3, :ntype2, :ntype1, :ntype2, :ntype3, :ntype1, :ntype3]
#     config = (solver=HydroModels.ODESolver(saveat=timeidx), ptypes=node_types, stypes=node_types, timeidx=timeidx)
#     output_arr = route(input_arr, pas, config=config)
#     return output_arr
# end


