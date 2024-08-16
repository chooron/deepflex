@testset "test simple flux (build by Symbolics.jl)" begin
    @variables a b c
    @parameters p1 p2
    simple_flux_1 = LumpedHydro.SimpleFlux([a, b] => [c], [p1, p2], exprs=[a * p1 + b * p2])
    @test LumpedHydro.get_input_names(simple_flux_1) == [:a, :b]
    @test LumpedHydro.get_param_names(simple_flux_1) == [:p1, :p2]
    @test LumpedHydro.get_output_names(simple_flux_1) == [:c,]
    @test simple_flux_1([2.0, 3.0], [3.0, 4.0]) == [2.0 * 3.0 + 3.0 * 4.0]
end

@testset "test simple flux (build by symbol)" begin
    simple_flux_1 = LumpedHydro.SimpleFlux([:a, :b] => [:c, :d], [:p1, :p2],
        flux_funcs=[(i, p) -> i[1] * p[1] + i[2] * p[2], (i, p) -> i[1] / p[1] + i[2] / p[2]])
    @test LumpedHydro.get_input_names(simple_flux_1) == [:a, :b]
    @test LumpedHydro.get_param_names(simple_flux_1) == [:p1, :p2]
    @test LumpedHydro.get_output_names(simple_flux_1) == [:c, :d]
    @test simple_flux_1([2.0, 3.0], [3.0, 4.0]) == [2.0 * 3.0 + 3.0 * 4.0, 2.0 / 3.0 + 3.0 / 4.0]
end

@testset "test state flux (build by Symbolics.jl)" begin
    @variables a b c d e
    @parameters p1 p2
    simple_flux_1 = LumpedHydro.SimpleFlux([a, b] => [c], [p1, p2], exprs=[a * p1 + b * p2])
    simple_flux_2 = LumpedHydro.SimpleFlux([a, b] => [d], [p1, p2], exprs=[a / p1 + b / p2])
    state_flux_1 = LumpedHydro.StateFlux([a, b] => [c, d], e)
    @test LumpedHydro.get_input_names(state_flux_1) == [:a, :b, :c, :d]
    @test LumpedHydro.get_param_names(state_flux_1) == Symbol[]
    @test LumpedHydro.get_output_names(state_flux_1) == Symbol[]
    @test LumpedHydro.get_state_names(state_flux_1) == [:e]

    state_flux_2 = LumpedHydro.StateFlux([a, b, c, d], e, [p1, p2], expr=a * p1 + b * p2 - c - d)
    @test LumpedHydro.get_input_names(state_flux_2) == [:a, :b, :c, :d]
    @test LumpedHydro.get_param_names(state_flux_2) == [:p1, :p2]
    @test LumpedHydro.get_output_names(state_flux_2) == Symbol[]
    @test LumpedHydro.get_state_names(state_flux_2) == [:e]

    state_flux_3 = LumpedHydro.StateFlux(d => e)
    @test LumpedHydro.get_input_names(state_flux_3) == [:e,]
    @test LumpedHydro.get_param_names(state_flux_3) == Symbol[]
    @test LumpedHydro.get_output_names(state_flux_3) == Symbol[]
    @test LumpedHydro.get_state_names(state_flux_3) == [:d,]
end

@testset "test neural flux" begin
    @variables a b c d e
    nn = Lux.Chain(
        layer_1=Lux.Dense(3, 16, Lux.leakyrelu),
        layer_2=Lux.Dense(16, 16, Lux.leakyrelu),
        layer_3=Lux.Dense(16, 2, Lux.leakyrelu),
        name=:testnn
    )
    nn_ps = LuxCore.initialparameters(StableRNG(42), nn)
    func = (x, p) -> LuxCore.stateless_apply(nn, x, p)
    nn_flux_1 = LumpedHydro.NeuralFlux([a, b, c] => [d, e], nn)
    @test LumpedHydro.get_input_names(nn_flux_1) == [:a, :b, :c]
    @test LumpedHydro.get_param_names(nn_flux_1) == Symbol[]
    @test LumpedHydro.get_nn_names(nn_flux_1) == [:testnn]
    @test LumpedHydro.get_output_names(nn_flux_1) == [:d, :e]
    @test nn_flux_1([1, 2, 3], Vector(ComponentVector(nn_ps))) ≈ func([1, 2, 3], nn_ps)
end