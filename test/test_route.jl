include("../src/HydroModels.jl")    

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

ndtypes = [:ntype1, :ntype2, :ntype3]
rflux = HydroModels.DischargeRouteFlux(q1)

params = ComponentVector(NamedTuple{Tuple(ndtypes)}([(lag=0.2,) for _ in eachindex(ndtypes)]))
initstates = ComponentVector(NamedTuple{Tuple(ndtypes)}([(s_river=0.1,) for _ in eachindex(ndtypes)]))
pas = ComponentVector(; params, initstates)
route = HydroModels.VectorRoute(:vectorroute; rfunc=rflux, network=network)

input_arr = ones(1, 9, 20)
timeidx = collect(1:20)
node_types = [:ntype1, :ntype2, :ntype3, :ntype2, :ntype1, :ntype2, :ntype3, :ntype1, :ntype3]
output_arr = route(input_arr, pas, timeidx=timeidx, ptypes=node_types)
# [params[ntype][:k] for ntype in node_types] .* input_arr[1,:,1]