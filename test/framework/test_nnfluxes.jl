using Lux
using LuxCore
using Symbolics
using StableRNGs
using ModelingToolkit
using ModelingToolkitNeuralNets
using ModelingToolkit: t_nounits as t
using ComponentArrays
using ModelingToolkitStandardLibrary.Blocks: RealInputArray
using RuntimeGeneratedFunctions
using BenchmarkTools

nn = Lux.Chain(
    Lux.Dense(4 => 16), # , Lux.tanh
    Lux.Dense(16 => 16), # , Lux.leakyrelu
    Lux.Dense(16 => 2)#, Lux.leakyrelu
)

init_params = Lux.initialparameters(StableRNG(42), nn)
init_states = Lux.initialstates(StableRNG(42), nn)
@parameters p[1:length(init_params)] = Vector(ComponentVector(init_params))
@parameters ptype::typeof(typeof(init_params)) = typeof(init_params) [tunable = false]
lazyconvert_p = Symbolics.array_term(convert, ptype, p, size=size(p))

@variables v1(t), v2(t), v3(t), v4(t)
@variables v(t)[1:4]
exprs = LuxCore.stateless_apply(nn, v, lazyconvert_p)[1]

temp_func = build_function(exprs, [v], [p, ptype], expression=Val{false})

temp_func([[3, 2, 1, 2]], [ComponentVector(init_params), typeof(ComponentVector(init_params))])


temp_func = build_function(v, [v], expression=Val{true})[2]

f = eval(temp_func)
f([[1,2,2,3]])


expr1 = v2 * 2
expr2 = v3 + 1

