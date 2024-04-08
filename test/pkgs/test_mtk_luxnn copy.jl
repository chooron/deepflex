using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Optimization, OptimizationOptimisers
using Lux
using OrdinaryDiffEq
using Random

model = Lux.Chain(Lux.Dense(3, 16, tanh), Lux.Dense(16, 1, leakyrelu), name=:model)
rng = MersenneTwister()
Random.seed!(rng, 42)
ps, st = Lux.setup(rng, model)

@variables x(t)
@parameters layer1_w[1:16, 1:3], layer1_b[1:16, 1:1]
@parameters layer2_w[1:16, 1:16], layer2_b[1:16, 1:1]

temp_func(t) = (model([sin(t), cos(t), tan(t)], (layer_1=(weight=layer1_w, bias=layer1_b), layer_2=(weight=layer2_w, bias=layer2_b)), st)[1])[1]

eqs = [D(x) ~ temp_func(t) .* x]
u0 = [x => 1.0]
pset = [
    layer1_w => ps[:layer_1][:weight],
    layer1_b => ps[:layer_1][:bias],
    layer2_w => ps[:layer_2][:weight],
    layer2_b => ps[:layer_2][:bias],
]
sys = ODESystem(eqs, t, name=:sys)
simple_sys = structural_simplify(sys)
prob = ODEProblem(simple_sys, u0, (0.0, 10.0), pset)
sol = solve(prob, Tsit5())



p, re = destructure(ps)
Lux.apply(model,[1,1,1],re(p),st)
# @parameters pvec[1:Lux.parameterlength(model)]
# @parameters pvec[1:Lux.parameterlength(model)]
# temp_func(t) = (model([sin(t), cos(t), tan(t)], re(pvec), st)[1])[1]
# eqs = [D(x) ~ temp_func(t) * (x)]
# pvec => p,
