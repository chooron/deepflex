using ModelingToolkit
using Graphs
using ComponentArrays
include("../src/HydroModels.jl")

@variables q1 q1_routed s_river q_gen
@parameters lag
rflux = HydroModels.HydroFlux([q1, s_river] => [q_gen, q1_routed], [lag], exprs=[q1, s_river / (1 + lag) + q1])

# network = DiGraph(9)
# add_edge!(network, 1, 2)
# add_edge!(network, 2, 5)
# add_edge!(network, 3, 5)
# add_edge!(network, 4, 5)
# add_edge!(network, 5, 8)
# add_edge!(network, 6, 9)
# add_edge!(network, 7, 8)
# add_edge!(network, 8, 9)

# ndtypes = [:ntype1, :ntype2, :ntype3]
# hrunames = [:nid1, :nid2, :nid3, :nid4, :nid5, :nid6, :nid7, :nid8, :nid9]
# params = ComponentVector(lag=fill(0.2, length(ndtypes)))
# initstates = ComponentVector(s_river=fill(0.0, length(hrunames)))
# pas = ComponentVector(; params, initstates)

# ptypes = [:ntype1, :ntype2, :ntype3, :ntype2, :ntype1, :ntype2, :ntype3, :ntype1, :ntype3]
# vroute = HydroModels.VectorRoute(rfunc=rflux, rstates=[s_river], network=network)

# ptyidx = [findfirst(isequal(ptype), ndtypes) for ptype in ptypes]
# styidx = [findfirst(isequal(hruname), hrunames) for hruname in hrunames]

# # 24 * 3600 / (10.0 * 1e6) * 1e3
# input_arr = rand(1, 9, 20)
# timeidx = collect(1:20)
# sol_2 = vroute(input_arr, pas, config=(timeidx=timeidx, ptyidx=ptyidx, styidx=styidx, solver=HydroModels.ManualSolver{true}()))


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
params = ComponentVector(rapid_k=fill(2.0, length(ndtypes)), rapid_x=fill(0.2, length(ndtypes)))
pas = ComponentVector(; params)
vroute = HydroModels.RapidRoute([q1]=>[q1_routed], network=network)
input_arr = ones(1, 9, 20)
timeidx = collect(1:20)
ptypes = [:ntype1, :ntype2, :ntype3, :ntype2, :ntype1, :ntype2, :ntype3, :ntype1, :ntype3]
sol_2 = vroute(input_arr, pas, config=(timeidx=timeidx, ptyidx=ptyidx))