using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Optimization, OptimizationOptimisers
using Flux
using OrdinaryDiffEq
using Random

model = Flux.Chain(Flux.Dense(3 => 16, Flux.tanh), Flux.Dense(16 => 1, Flux.leakyrelu))
# rng = MersenneTwister()
# Random.seed!(rng, 42)
# ps, st = Lux.setup(rng, model)

ps, re = destructure(model)
re(ps)([1, 1, 1])

@variables x(t)
@parameters pvec[1:Lux.parameterlength(model)]

temp_func(t) = re(pvec)([sin(t), cos(t), tan(t)])[1]
eqs = [D(x) ~ temp_func(t) * (x)]
u0 = [x => 1.0]
pset = [
    pvec => ps,
]
sys = ODESystem(eqs, t, name=:sys)
simple_sys = structural_simplify(sys)
prob = ODEProblem(simple_sys, u0, (0.0, 10.0), pset)
sol = solve(prob, Tsit5())

