include("../src/HydroModels.jl")
using ModelingToolkit
using ComponentArrays
# @variables q1

# flwdir = [1 4 8; 1 4 4; 1 1 2]
# positions = [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]

# ndtypes = [:ntype1, :ntype2, :ntype3]
# rflux = HydroModels.CascadeRouteFlux(q1)

# params = ComponentVector(NamedTuple{Tuple(ndtypes)}([(n=3, k=0.2) for _ in eachindex(ndtypes)]))
# pas = ComponentVector(; params)
# route = HydroModels.GridRoute(:gridroute; rfunc=rflux, flwdir=flwdir, positions=positions)

# input_arr = ones(1, 9, 20)
# timeidx = collect(1:20)
# node_types = [:ntype1, :ntype2, :ntype3, :ntype2, :ntype1, :ntype2, :ntype3, :ntype1, :ntype3]
# output_arr = route(input_arr, pas, timeidx, node_types)

@variables q1

flwdir = [1 4 8; 1 4 4; 1 1 2]
positions = [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]

ndtypes = [:ntype1, :ntype2, :ntype3]
rflux = HydroModels.CascadeRouteFlux(q1)

params = ComponentVector(NamedTuple{Tuple(ndtypes)}([(n=3, k=0.2) for _ in eachindex(ndtypes)]))
pas = ComponentVector(; params)
route = HydroModels.VectorRoute(:gridroute; rfunc=rflux, flwdir=flwdir, positions=positions)

input_arr = ones(1, 9, 20)
timeidx = collect(1:20)
node_types = [:ntype1, :ntype2, :ntype3, :ntype2, :ntype1, :ntype2, :ntype3, :ntype1, :ntype3]
output_arr = route(input_arr, pas, timeidx, node_types)