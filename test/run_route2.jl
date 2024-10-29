include("../src/HydroModels.jl")
using ModelingToolkit
using ComponentArrays
using Graphs
using Plots

solver_1 = HydroModels.DiscreteSolver()
solver_2 = HydroModels.ODESolver()


function run_vector_route()
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
    params = ComponentVector(NamedTuple{Tuple(ndtypes)}([(k=2.0, x=0.2) for _ in eachindex(ndtypes)]))
    initstates = ComponentVector(NamedTuple{Tuple(ndtypes)}([(s_river=0.0,) for _ in eachindex(ndtypes)]))
    pas = ComponentVector(; params, initstates)
    vroute = HydroModels.VectorRoute(rfunc=rflux_1, network=network, subareas=10.0)
    # 24 * 3600 / (10.0 * 1e6) * 1e3
    input_arr = ones(1, 9, 20)
    timeidx = collect(1:20)
    ptypes = [:ntype1, :ntype2, :ntype3, :ntype2, :ntype1, :ntype2, :ntype3, :ntype1, :ntype3]
    sol_2 = vroute(input_arr, pas, config=(timeidx=timeidx, ptypes=ptypes, solver=solver_1))
    return sol_2
end

function run_grid_route()
    @variables q
    rflux = HydroModels.DischargeRouteFlux(q)
    flwdir = [1 4 8; 1 4 4; 1 1 2]
    positions = [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]

    params = ComponentVector(NamedTuple{Tuple(ndtypes)}([(lag=0.2,) for _ in eachindex(ndtypes)]))
    initstates = ComponentVector(NamedTuple{Tuple(ndtypes)}([(s_river=0.1,) for _ in eachindex(ndtypes)]))
    pas = ComponentVector(; params, initstates)
    route = HydroModels.GridRoute(rfunc=rflux, flwdir=flwdir, positions=positions, subareas=10.0)

    input_arr = ones(1, 9, 20)
    timeidx = collect(1:20)
    node_types = [:ntype1, :ntype2, :ntype3, :ntype2, :ntype1, :ntype2, :ntype3, :ntype1, :ntype3]
    config = (solver=HydroModels.ODESolver(saveat=timeidx), ptypes=node_types, stypes=node_types, timeidx=timeidx)
    output_arr = route(input_arr, pas, config=config)
    return output_arr
end

output_arr_1 = run_grid_route()
plot(output_arr_1[2, :, :]')
