include("../src/HydroModels.jl")    
using ModelingToolkit
using Graphs
using ComponentArrays
@variables q1

network = DiGraph(9)
add_edge!(network, 1, 2)
add_edge!(network, 2, 5)
add_edge!(network, 3, 5)
add_edge!(network, 4, 5)
add_edge!(network, 5, 8)
add_edge!(network, 6, 9)
add_edge!(network, 7, 8)
add_edge!(network, 8, 9)
nv(network)
ndtypes = [:ntype1, :ntype2, :ntype3]
rflux = HydroModels.DischargeRouteFlux(q1)

params = ComponentVector(NamedTuple{Tuple(ndtypes)}([(lag=0.2,) for _ in eachindex(ndtypes)]))
initstates = ComponentVector(NamedTuple{Tuple(ndtypes)}([(s_river=0.1,) for _ in eachindex(ndtypes)]))
pas = ComponentVector(; params, initstates)
route = HydroModels.VectorRoute(name=:vectorroute, rfunc=rflux, network=network, subareas=10.0)

input_arr = ones(1, 9, 20)
timeidx = (collect(1:20) .- 1 ).*24
node_types = [:ntype1, :ntype2, :ntype3, :ntype2, :ntype1, :ntype2, :ntype3, :ntype1, :ntype3]
output_arr = route(input_arr, pas, timeidx, ptypes=node_types, config=(solver=HydroModels.ODESolver(saveat=timeidx), interp=LinearInterpolation, ptypes=node_types))
# [params[ntype][:k] for ntype in node_types] .* input_arr[1,:,1]