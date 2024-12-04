@testset "test unit hydro flux" begin
    # Define the variables and parameters
    @variables q1 q1_routed
    @parameters x1
    # Define a unit hydrograph function
    uh_func = HydroModels.UHFunction(:UH_1_HALF)
    # Create a UnitHydroRouteFlux object
    # Input: q1 (flow)
    # Parameter: x1 (routing parameter)
    # Using uh_1_half as the unit hydrograph function
    # Solve type: unithydro1 (convolution method)
    uh1 = HydroModels.UnitHydrograph(q1, q1_routed, x1, uhfunc=uh_func, solvetype=:DISCRETE)
    uh2 = HydroModels.UnitHydrograph(q1, q1_routed, x1, uhfunc=uh_func, solvetype=:SPARSE)
    # Test the input names of the router
    @test HydroModels.get_input_names(uh1) == HydroModels.get_input_names(uh2) == [:q1]
    # Test the parameter names of the router
    @test HydroModels.get_param_names(uh1) == HydroModels.get_param_names(uh2) == [:x1]
    # Test the output names of the router
    @test HydroModels.get_output_names(uh1) == HydroModels.get_output_names(uh2) == [:q1_routed]
    # Test the routing function with sample input
    input_flow = Float32[2 3 4 2 3 1]
    params = ComponentVector(params=(x1=3.5,))
    expected_output = [0.0899066  0.643448  2.3442  3.20934  3.44646  2.20934]
    # [0.08726897695099571 0.5373023715895905 1.6508571480656808 2.839759323622619 3.2301609643779736 2.7991762465729138]
    @test uh1(input_flow, params) ≈ expected_output atol = 1e-3
    @test uh2(input_flow, params) ≈ expected_output atol = 1e-3
    # test with multiple nodes
    input_arr = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), input_flow for _ in 1:10), (1, 3, 2))
    ndtypes = [Symbol("node_$i") for i in 1:10]
    node_params = NamedTuple{Tuple(ndtypes)}(repeat([(x1=3.5,)], 10))
    node_initstates = NamedTuple{Tuple(ndtypes)}(repeat([NamedTuple()], 10))
    input_pas = ComponentVector(params=node_params, initstates=node_initstates)

    @test uh1(input_arr, input_pas, ptypes=ndtypes) ≈ permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), expected_output for _ in 1:10), (1, 3, 2)) atol = 1e-3
    @test uh2(input_arr, input_pas, ptypes=ndtypes) ≈ permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), expected_output for _ in 1:10), (1, 3, 2)) atol = 1e-3
end