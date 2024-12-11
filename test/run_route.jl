@testset "test grid route based on discharge route flux" begin
    @variables q1 q1_routed s_river
    @parameters lag
    flwdir = [1 4 8; 1 4 4; 1 1 2]
    positions = [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]
    
    ndtypes = [:ntype1, :ntype2, :ntype3]
    hrunames = [:nid1, :nid2, :nid3, :nid4, :nid5, :nid6, :nid7, :nid8, :nid9]
    rflux = HydroModels.HydroFlux([q1, s_river] => [q1_routed], [lag], exprs=[s_river / (1 + lag) + q1])
    @test HydroModels.get_input_names(rflux) == [:q1, :s_river]
    @test HydroModels.get_output_names(rflux) == [:q1_routed]
    @test HydroModels.get_param_names(rflux) == [:lag]
    params = ComponentVector(NamedTuple{Tuple(ndtypes)}([(lag=0.2,) for _ in eachindex(ndtypes)]))
    initstates = ComponentVector(NamedTuple{Tuple(hrunames)}([(s_river=0.1,) for _ in eachindex(hrunames)]))
    pas = ComponentVector(; params, initstates)
    route = HydroModels.GridRoute(rfunc=rflux, rstate=s_river, flwdir=flwdir, positions=positions)
    @test HydroModels.get_input_names(route) == [:q1]
    @test HydroModels.get_output_names(route) == [:q1_routed]
    @test HydroModels.get_param_names(route) == [:lag]
    @test HydroModels.get_state_names(route) == [:s_river]
    
    input_arr = ones(1, 9, 20)
    timeidx = collect(1:20)
    node_types = [:ntype1, :ntype2, :ntype3, :ntype2, :ntype1, :ntype2, :ntype3, :ntype1, :ntype3]
    config = (solver=ManualSolver{true}(), interp=LinearInterpolation, ptypes=node_types, stypes=hrunames, timeidx=timeidx)
    output_arr = route(input_arr, pas, config=config)
    @test size(output_arr) == size(ones(2, 9, 20))
end

@testset "test vector route based on discharge route flux" begin
    @variables q1 q1_routed s_river
    @parameters lag
    rflux = HydroModels.HydroFlux([q1, s_river] => [q1_routed], [lag], exprs=[s_river / (1 + lag) + q1])

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
    hrunames = [:nid1, :nid2, :nid3, :nid4, :nid5, :nid6, :nid7, :nid8, :nid9]
    params = ComponentVector(NamedTuple{Tuple(ndtypes)}([(lag=0.2,) for _ in eachindex(ndtypes)]))
    initstates = ComponentVector(NamedTuple{Tuple(hrunames)}([(s_river=0.0,) for _ in eachindex(hrunames)]))
    pas = ComponentVector(; params, initstates)
    vroute = HydroModels.VectorRoute(rfunc=rflux, rstate=s_river, network=network)
    # 24 * 3600 / (10.0 * 1e6) * 1e3
    input_arr = ones(1, 9, 20)
    timeidx = collect(1:20)
    ptypes = [:ntype1, :ntype2, :ntype3, :ntype2, :ntype1, :ntype2, :ntype3, :ntype1, :ntype3]
    sol_2 = vroute(input_arr, pas, config=(timeidx=timeidx, ptypes=ptypes, stypes=hrunames, solver=ManualSolver{true}()))
    @test size(sol_2) == size(ones(2, 9, 20))
end

@testset "test rapid route based on muskingum route flux" begin
    @variables q q_routed

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
    params = ComponentVector(NamedTuple{Tuple(ndtypes)}([(rapid_k=2.0, rapid_x=0.2) for _ in eachindex(ndtypes)]))
    initstates = ComponentVector(NamedTuple{Tuple(ndtypes)}([(s_river=0.0,) for _ in eachindex(ndtypes)]))
    pas = ComponentVector(; params, initstates)
    vroute = HydroModels.RapidRoute([q]=>[q_routed], network=network)
    input_arr = ones(1, 9, 20)
    timeidx = collect(1:20)
    ptypes = [:ntype1, :ntype2, :ntype3, :ntype2, :ntype1, :ntype2, :ntype3, :ntype1, :ntype3]
    sol_2 = vroute(input_arr, pas, config=(timeidx=timeidx, ptypes=ptypes, solver=ManualSolver{true}()))
    @test size(sol_2) == size(ones(1, 9, 20))
end
