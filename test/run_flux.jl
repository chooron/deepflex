@testset "test hydro flux (build by Symbolics.jl)" begin
    @variables a b c
    @parameters p1 p2
    simple_flux_1 = HydroModels.HydroFlux([a, b] => [c], [p1, p2], exprs=[a * p1 + b * p2])
    @test HydroModels.get_input_names(simple_flux_1) == [:a, :b]
    @test HydroModels.get_param_names(simple_flux_1) == [:p1, :p2]
    @test HydroModels.get_output_names(simple_flux_1) == [:c,]
    @test simple_flux_1([2.0, 3.0], ComponentVector(params=(p1=3.0, p2=4.0))) == [2.0 * 3.0 + 3.0 * 4.0]
    output_mat = [18.0 17.0 11.0]
    @test simple_flux_1([2.0 3.0 1.0; 3.0 2.0 2.0], ComponentVector(params=(p1=3.0, p2=4.0))) == output_mat
    # test with multiple nodes
    input_arr = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), [2.0 3.0 1.0; 3.0 2.0 2.0] for _ in 1:10), (1, 3, 2))
    ndtypes = [Symbol("node_$i") for i in 1:10]
    input_pas = ComponentVector(params=NamedTuple{Tuple(ndtypes)}(repeat([(p1=3.0, p2=4.0)], 10)))
    @test simple_flux_1(input_arr, input_pas, ptypes=ndtypes) == permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), output_mat for _ in 1:10), (3, 1, 2))
end

@testset "test state flux (build by Symbolics.jl)" begin
    # Define variables and parameters
    @variables a b c d e
    @parameters p1 p2

    # Create HydroFlux objects for comparison
    simple_flux_1 = HydroModels.HydroFlux([a, b] => [c], [p1, p2], exprs=[a * p1 + b * p2])
    simple_flux = HydroModels.HydroFlux([a, b] => [d], [p1, p2], exprs=[a / p1 + b / p2])

    # Test StateFlux with input and output, but no parameters
    state_flux_1 = HydroModels.StateFlux([a, b] => [c, d], e)
    @test HydroModels.get_input_names(state_flux_1) == [:a, :b, :c, :d]
    @test HydroModels.get_param_names(state_flux_1) == Symbol[]
    @test HydroModels.get_output_names(state_flux_1) == Symbol[]
    @test HydroModels.get_state_names(state_flux_1) == [:e]

    # Test StateFlux with input, parameters, and a custom expression
    state_flux = HydroModels.StateFlux([a, b, c, d], e, [p1, p2], expr=a * p1 + b * p2 - c - d)
    @test HydroModels.get_input_names(state_flux) == [:a, :b, :c, :d]
    @test HydroModels.get_param_names(state_flux) == [:p1, :p2]
    @test HydroModels.get_output_names(state_flux) == Symbol[]
    @test HydroModels.get_state_names(state_flux) == [:e]

    # Test StateFlux with a single input and state
    state_flux_3 = HydroModels.StateFlux(d => e)
    @test HydroModels.get_input_names(state_flux_3) == [:e,]
    @test HydroModels.get_param_names(state_flux_3) == Symbol[]
    @test HydroModels.get_output_names(state_flux_3) == Symbol[]
    @test HydroModels.get_state_names(state_flux_3) == [:d,]
end

# todo muskingum need rebuild
# @testset "test muskingum route flux" begin
#     @variables q1

#     # Building the Muskingum routing flux
#     k, x = 3.0, 0.2
#     pas = ComponentVector(params=(k=k, x=x,))
#     msk_flux = HydroModels.MuskingumRouteFlux(q1)
#     input = Float64[1 2 3 2 3 2 5 7 8 3 2 1]
#     re = msk_flux(input, pas)

#     # Verifying the input, output, and parameter names
#     @test HydroModels.get_input_names(msk_flux) == [:q1]
#     @test HydroModels.get_output_names(msk_flux) == [:q1_routed]
#     @test HydroModels.get_param_names(msk_flux) == [:k, :x]

#     # Checking the size and values of the output
#     @test size(re) == size(input)
#     @test re ≈ [1.0 0.977722 1.30086 1.90343 1.919 2.31884 2.15305 3.07904 4.39488 5.75286 4.83462 3.89097] atol = 1e-1
# end

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
    uhflux1 = HydroModels.UnitHydroFlux(q1, q1_routed, x1, uhfunc=uh_func, solvetype=:DISCRETE)
    uhflux2 = HydroModels.UnitHydroFlux(q1, q1_routed, x1, uhfunc=uh_func, solvetype=:SPARSE)
    # Test the input names of the router
    @test HydroModels.get_input_names(uhflux1) == HydroModels.get_input_names(uhflux2) == [:q1]
    # Test the parameter names of the router
    @test HydroModels.get_param_names(uhflux1) == HydroModels.get_param_names(uhflux2) == [:x1]
    # Test the output names of the router
    @test HydroModels.get_output_names(uhflux1) == HydroModels.get_output_names(uhflux2) == [:q1_routed]
    # Test the routing function with sample input
    input_flow = Float32[2 3 4 2 3 1]
    params = ComponentVector(params=(x1=3.5,))
    expected_output = [0.08990658541313483 0.6434483278713437 2.3442008102889558 3.209340932169254 3.4464582575417912 2.2093409321692534]
    @test uhflux1(input_flow, params) ≈ expected_output
    @test uhflux2(input_flow, params) ≈ expected_output
    # test with multiple nodes
    input_arr = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), input_flow for _ in 1:10), (1, 3, 2))
    ndtypes = [Symbol("node_$i") for i in 1:10]
    node_params = NamedTuple{Tuple(ndtypes)}(repeat([(x1=3.5,)], 10))
    node_initstates = NamedTuple{Tuple(ndtypes)}(repeat([NamedTuple()], 10))
    input_pas = ComponentVector(params=node_params, initstates=node_initstates)
    @test uhflux1(input_arr, input_pas, ptypes=ndtypes) ≈ permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), expected_output for _ in 1:10), (1, 3, 2))
    @test uhflux2(input_arr, input_pas, ptypes=ndtypes) ≈ permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), expected_output for _ in 1:10), (1, 3, 2))
end

@testset "test neural flux (single output)" begin
    @variables a b c d e
    nn_1 = Lux.Chain(
        layer_1=Lux.Dense(3, 16, Lux.leakyrelu),
        layer=Lux.Dense(16, 16, Lux.leakyrelu),
        layer_3=Lux.Dense(16, 1, Lux.leakyrelu),
        name=:testnn
    )
    nn_ps = ComponentVector(LuxCore.initialparameters(StableRNG(42), nn_1))
    nn_ps_vec = collect(ComponentVector(LuxCore.initialparameters(StableRNG(42), nn_1)))
    func_1 = (x, p) -> LuxCore.stateless_apply(nn_1, x, p)
    nn_flux_1 = HydroModels.NeuralFlux([a, b, c] => [d], nn_1)
    @test HydroModels.get_input_names(nn_flux_1) == [:a, :b, :c]
    @test HydroModels.get_param_names(nn_flux_1) == Symbol[]
    @test HydroModels.get_nn_names(nn_flux_1) == [:testnn]
    @test HydroModels.get_output_names(nn_flux_1) == [:d]
    @test nn_flux_1([1, 2, 3], ComponentVector(nn=(testnn=nn_ps_vec,))) ≈ func_1([1, 2, 3], nn_ps)
    test_input = [1 3 3; 2 2 2; 1 2 1; 3 1 2]
    expected_output = func_1(test_input', nn_ps)
    @test nn_flux_1(test_input, ComponentVector(nn=(testnn=nn_ps_vec,))) ≈ expected_output
    # test with multiple nodes
    input_arr = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), test_input for _ in 1:10), (2, 3, 1))
    input_pas = ComponentVector(nn=(testnn=nn_ps_vec,))
    @test nn_flux_1(input_arr, input_pas) ≈ permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), expected_output for _ in 1:10), (1, 3, 2))
end

@testset "test neural flux (multiple output)" begin
    @variables a b c d e
    nn = Lux.Chain(
        layer_1=Lux.Dense(3, 16, Lux.leakyrelu),
        layer_2=Lux.Dense(16, 16, Lux.leakyrelu),
        layer_3=Lux.Dense(16, 2, Lux.leakyrelu),
        name=:testnn
    )
    nn_ps = LuxCore.initialparameters(StableRNG(42), nn)
    nn_ps_vec = collect(ComponentVector(nn_ps))
    func = (x, p) -> LuxCore.stateless_apply(nn, x, p)
    nn_flux = HydroModels.NeuralFlux([a, b, c] => [d, e], nn)
    @test HydroModels.get_input_names(nn_flux) == [:a, :b, :c]
    @test HydroModels.get_param_names(nn_flux) == Symbol[]
    @test HydroModels.get_nn_names(nn_flux) == [:testnn]
    @test HydroModels.get_output_names(nn_flux) == [:d, :e]
    test_input = [1 3 3; 2 2 2; 1 2 1; 3 1 2]
    test_output = func(test_input', nn_ps)
    @test nn_flux([1, 2, 3], ComponentVector(nn=(testnn=nn_ps_vec,))) ≈ func([1, 2, 3], nn_ps)
    @test nn_flux(test_input, ComponentVector(nn=(testnn=nn_ps_vec,))) ≈ test_output
    # test with multiple nodes
    input_arr = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), test_input for _ in 1:10), (2, 3, 1))
    input_pas = ComponentVector(nn=(testnn=nn_ps_vec,))
    @test nn_flux(input_arr, input_pas) ≈ permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), test_output for _ in 1:10), (1, 3, 2))
end

@testset "test runtime flux function building" begin
    @variables a b c d e
    @parameters x1 x2
    inputs = [a, b, c]
    outputs = [d, e]
    params = [x1, x2]
    exprs = [x1 * a + x2 * b + c, x1 * c + x2]
    func = HydroModels.build_flux_func(inputs, outputs, params, exprs)
    @test func([1, 2, 1], [2, 1], 1) == [2 * 1 + 1 * 2 + 1, 2 * 1 + 1]
end