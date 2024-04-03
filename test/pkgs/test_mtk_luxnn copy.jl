using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Optimization, OptimizationOptimisers
using Lux
using OrdinaryDiffEq
using Random

model = Lux.Chain(Lux.Dense(3, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 1, leakyrelu), name=:model)
rng = MersenneTwister()
Random.seed!(rng, 42)
ps, st = Lux.setup(rng, model)
time = 1:100
x1 = sin.(time)
x2 = cos.(time)
x3 = tan.(time)
x = hcat(x1, x2, x3)'
y = @.(2 * x1 + 3 * x2 - 5 * x3) + rand(100)
y_hat = vec(model([1, 2, 3], ps, st)[1])

@variables x(t)
@parameters layer1_w[1:16, 1:3] layer1_b[1:16, 1:1]
@parameters layer2_w[1:16, 1:16] layer2_b[1:16, 1:1]
@parameters layer3_w[1:1, 1:16] layer3_b[1:1, 1:1]

p_tuple = (layer_1=(weight=layer1_w, bias=layer1_b), layer_2=(weight=layer2_w, bias=layer2_b), layer_3=(weight=layer3_w, bias=layer3_b))
temp_func(t) = (model([sin(t), cos(t), tan(t)], p_tuple, st)[1])[1]
eqs = [D(x) ~ temp_func(t) / (x)]
u0 = [x => 1.0]
p = [
    layer1_w => ps[:layer_1][:weight],
    layer1_b => ps[:layer_1][:bias],
    layer2_w => ps[:layer_2][:weight],
    layer2_b => ps[:layer_2][:bias],
    layer3_w => ps[:layer_3][:weight],
    layer3_b => ps[:layer_3][:bias],
]
sys = ODESystem(eqs, t, name=:sys)
simple_sys = structural_simplify(sys)
prob = ODEProblem(simple_sys, u0, (0.0, 10.0), p)
sol = solve(prob,Tsit5())

