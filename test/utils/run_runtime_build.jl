@testset "test runtime flux function building" begin
    @variables a b c d e
    @parameters x1 x2
    inputs = [a, b, c]
    outputs = [d, e]
    params = [x1, x2]
    exprs = [x1 * a + x2 * b + c, x1 * c + x2]
    func = HydroModels.build_flux_func(inputs, outputs, params, exprs)
    @test func([1, 2, 1], [2, 1]) == [2 * 1 + 1 * 2 + 1, 2 * 1 + 1]
end

@testset "test build state function" begin
    @variables routingstore exch slowflow_routed routedflow
    @parameters x2, x3
    funcs = [
        SimpleFlux([routingstore] => [exch], [x2, x3],
            exprs=[x2 * abs(routingstore / x3)^3.5]),
        SimpleFlux([routingstore, slowflow_routed, exch] => [routedflow], [x3],
            exprs=[x3^(-4) / 4 * (routingstore + slowflow_routed + exch)^5]),
    ]
    dfunc = HydroModels.StateFlux([slowflow_routed, exch] => [routedflow], routingstore)
    infos = (input=[:slowflow_routed], state=[:routingstore])
    flux_func, state_func = HydroModels.build_ele_func(funcs, [dfunc], infos)
    rgt, slg = 10.0, 20.0
    exch = 2.42 * abs(rgt / 69.63)^3.5
    routedflow = 69.63^(-4) / 4 * (rgt + slg + exch)^5
    @test flux_func([slg,rgt],[2.42, 69.63], []) ≈  [exch, routedflow]
    @test state_func([slg], [rgt], [2.42, 69.63], []) ≈ [exch + slg - 69.63^(-4) / 4 * (rgt + slg + exch)^5]
end