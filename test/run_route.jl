@testset "test grid route based on hydrology discharge route flux" begin
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
    route = HydroModels.GridRoute(:gridroute; rfunc=rflux, flwdir=flwdir, positions=positions)
    @test HydroModels.get_input_names(route) == [:q1]
    @test HydroModels.get_output_names(route) == [:q1_routed]
    @test HydroModels.get_param_names(route) == [:lag]

    input_arr = ones(1, 9, 20)
    timeidx = collect(1:20)
    node_types = [:ntype1, :ntype2, :ntype3, :ntype2, :ntype1, :ntype2, :ntype3, :ntype1, :ntype3]
    output_arr = route(input_arr, pas, timeidx, node_types)
    #* we cannot get test data for now, thus we just test it is run success

    @test size(output_arr) == size(input_arr)
end

@testset "test grid route based on hydrology discharge route flux" begin
    @variables q1

    flwdir = [1 4 8; 1 4 4; 1 1 2]
    positions = [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]
    
    ndtypes = [:ntype1, :ntype2, :ntype3]
    rflux = HydroModels.CascadeRouteFlux(q1)
    
    params = ComponentVector(NamedTuple{Tuple(ndtypes)}([(n=3, k=0.2) for _ in eachindex(ndtypes)]))
    pas = ComponentVector(; params)
    route = HydroModels.GridRoute(:gridroute; rfunc=rflux, flwdir=flwdir, positions=positions)
    
    input_arr = ones(1, 9, 20)
    timeidx = collect(1:20)
    node_types = [:ntype1, :ntype2, :ntype3, :ntype2, :ntype1, :ntype2, :ntype3, :ntype1, :ntype3]
    output_arr = route(input_arr, pas, timeidx, node_types)

    @test size(output_arr) == size(input_arr)
end