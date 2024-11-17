@testset "test grid route based on discharge route flux" begin
    @variables q1 q1_routed s_river
    @parameters lag
    flwdir = [1 4 8; 1 4 4; 1 1 2]
    positions = [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]
    
    ndtypes = [:ntype1, :ntype2, :ntype3]
    rflux = HydroModels.HydroFlux([q1, s_river] => [q1_routed], [lag], exprs=[s_river / (1 + lag) + q1])
    @test HydroModels.get_input_names(rflux) == [:q1, :s_river]
    @test HydroModels.get_output_names(rflux) == [:q1_routed]
    @test HydroModels.get_param_names(rflux) == [:lag]
    params = ComponentVector(NamedTuple{Tuple(ndtypes)}([(lag=0.2,) for _ in eachindex(ndtypes)]))
    initstates = ComponentVector(NamedTuple{Tuple(ndtypes)}([(s_river=0.1,) for _ in eachindex(ndtypes)]))
    pas = ComponentVector(; params, initstates)
    nodeids = [:nid1, :nid2, :nid3, :nid4, :nid5, :nid6, :nid7, :nid8, :nid9]
    route = HydroModels.GridRoute(rfuncs=[rflux], rstates=[s_river], flwdir=flwdir, positions=positions, subareas=10.0, nodeids=nodeids)
    @test HydroModels.get_input_names(route) == [:q1]
    @test HydroModels.get_output_names(route) == [:q1_routed]
    @test HydroModels.get_param_names(route) == [:lag]
    @test HydroModels.get_state_names(route) == [:s_river]
    
    input_arr = ones(1, 9, 20)
    timeidx = collect(1:20)
    node_types = [:ntype1, :ntype2, :ntype3, :ntype2, :ntype1, :ntype2, :ntype3, :ntype1, :ntype3]
    config = (solver=HydroModels.ODESolver(saveat=timeidx), interp=LinearInterpolation, ptypes=node_types, stypes=node_types, timeidx=timeidx)
    output_arr = route(input_arr, pas, config=config)
    @test size(output_arr) == size(ones(2, 9, 20))

    @variables q1

    flwdir = [1 4 8; 1 4 4; 1 1 2]
    positions = [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]

    ndtypes = [:ntype1, :ntype2, :ntype3]
    rflux = HydroModels.DischargeRouteFlux(q1)
    @test HydroModels.get_input_names(rflux) == [:q1]
    @test HydroModels.get_output_names(rflux) == [:q1_routed]
    @test HydroModels.get_param_names(rflux) == [:lag]
    @test HydroModels.get_state_names(rflux) == [:s_river]
    params = ComponentVector(NamedTuple{Tuple(ndtypes)}([(lag=0.2,) for _ in eachindex(ndtypes)]))
    initstates = ComponentVector(NamedTuple{Tuple(ndtypes)}([(s_river=0.1,) for _ in eachindex(ndtypes)]))
    pas = ComponentVector(; params, initstates)
    route = HydroModels.GridRoute(rfunc=rflux, flwdir=flwdir, positions=positions, subareas=10.0)
    @test HydroModels.get_input_names(route) == [:q1]
    @test HydroModels.get_output_names(route) == [:q1_routed]
    @test HydroModels.get_param_names(route) == [:lag]

    input_arr = ones(1, 9, 20)
    timeidx = collect(1:20)
    node_types = [:ntype1, :ntype2, :ntype3, :ntype2, :ntype1, :ntype2, :ntype3, :ntype1, :ntype3]
    config = (solver=HydroModels.ODESolver(saveat=timeidx), interp=LinearInterpolation, ptypes=node_types, stypes=node_types, timeidx=timeidx)
    output_arr = route(input_arr, pas, config=config)
    #* we cannot get test data for now, thus we just test it is run success

    @test size(output_arr) == size(ones(2, 9, 20))
end

@testset "test vector route based on muskingum route flux" begin
    @variables q
    rflux_1 = HydroModels.RiverRouteFlux(q)

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
    sol_2 = vroute(input_arr, pas, config=(timeidx=timeidx, ptypes=ptypes, stypes=ptypes, solver=HydroModels.DiscreteSolver()))
    @test size(sol_2) == size(ones(2, 9, 20))
end


@testset "test rapid route based on muskingum route flux" begin
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
    vroute = HydroModels.RapidRoute(rfunc=rflux_1, network=network, subareas=10.0)
    # 24 * 3600 / (10.0 * 1e6) * 1e3
    input_arr = ones(1, 9, 20)
    timeidx = collect(1:20)
    ptypes = [:ntype1, :ntype2, :ntype3, :ntype2, :ntype1, :ntype2, :ntype3, :ntype1, :ntype3]
    sol_2 = vroute(input_arr, pas, config=(timeidx=timeidx, ptypes=ptypes, solver=HydroModels.DiscreteSolver()))
    @test size(sol_2) == size(ones(1, 9, 20))
end
