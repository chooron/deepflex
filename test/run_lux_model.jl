using Lux
using LuxCore
using Symbolics
using ModelingToolkit
using BenchmarkTools
model = Lux.Chain(
    Lux.Dense(3, 64),
    Lux.relu,
    Lux.Dense(64, 1),
    Lux.softmax
)

@variables a b c d
@parameters p1 p2 p3
model([a, b, c]=>[d])
(chain::LuxCore.AbstractExplicitContainerLayer)(var::Vector{Num}) = (chain=chain, input=var)
ModelingToolkit.isparameter(p1)
ModelingToolkit.isvariable(p1)
eq1 = a + b + c*p1 ~ d
vars = Num.(get_variables(eq1))
filter(x->!ModelingToolkit.isparameter(x), vars)
filter(x->ModelingToolkit.isparameter(x), vars)
