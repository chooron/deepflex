using ModelingToolkit
using Optimization, OptimizationOptimisers
using Lux
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
y_hat = vec(model(x, ps, st)[1])

@variables layer1_w[1:16, 1:3], layer1_b[1:16, 1:1]
@variables layer2_w[1:16, 1:16], layer2_b[1:16, 1:1]
@variables layer3_w[1:1, 1:16], layer3_b[1:1, 1:1]

# function loss_func(p)
#     y_hat = vec(model(x, p, st)[1])
#     sum(abs(y .- y_hat))
# end
p_tuple = (layer_1=(weight=layer1_w, bias=layer1_b), layer_2=(weight=layer2_w, bias=layer2_b), layer_3=(weight=layer3_w, bias=layer3_b))
loss_func = (sum((y .- vec(model(x, p_tuple, st)[1])) .^ 2))


sys = OptimizationSystem(loss_func, [layer1_w; layer1_b; layer2_w; layer2_b; layer3_w; layer3_b], [], name=:sys)
sys = complete(sys)
u0 = [
    layer1_w => ps[:layer_1][:weight],
    layer1_b => ps[:layer_1][:bias],
    layer2_w => ps[:layer_2][:weight],
    layer2_b => ps[:layer_2][:bias],
    layer3_w => ps[:layer_3][:weight],
    layer3_b => ps[:layer_3][:bias],
]
prob = OptimizationProblem(sys, u0)
solve(prob, GradientDescent())

