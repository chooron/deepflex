include("../src/HydroModels.jl")

using Symbolics

@variables q2 i2 q1 i1 g
@parameters c0 c1 c2 A

func1(x) = ones(3, 3) * x
@register_symbolic func1(x)

expr = c2 * q1 + c1 * i2 + c0 * i1
flux = HydroModels.HydroFlux([q1, i1, i2] => [q2], [c0, c1, c2], exprs=[expr])
flux1 = HydroModels.HydroFlux([q1, g] => [i1], exprs=[func1(q1 + g)])
flux2 = HydroModels.HydroFlux([q2, g] => [i2], exprs=[func1(q2 + g)])

flux2.func([[1, 2, 3], [2, 2, 1]], [], 1)


get_variables(func1(x))