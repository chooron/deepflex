using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)
using Symbolics
using SymbolicUtils
using SymbolicUtils.Code
using Lux
using LuxCore
using StableRNGs
using ModelingToolkit
using ComponentArrays

@variables a b c d e f(e)
@variables v[1:1]
@variables v1 v2 v3 v4
ass = [
    Assignment(a, b + c),
    Assignment(e, c + d),
    Assignment(v, MakeArray([a, b, c, d], Vector)),
    Assignment(e, v[1])
]

let_ = Let(ass, a + e + sum(v), false)
func = @RuntimeGeneratedFunction(
    toexpr(Func([a, e], [Assignment(var, 0.0) for var in [a, b, c, d, e]], let_))
)

func(1, 2)

chain = Lux.Chain(Lux.Dense(2, 16), Lux.Dense(16, 2))
init_params = ComponentVector(Lux.initialparameters(StableRNG(42), chain))
chain_params = first(@parameters p[1:length(init_params)] = Vector(ComponentVector(init_params)))
@parameters ptype::typeof(typeof(init_params)) = typeof(init_params) [tunable = false]
lazyconvert_params = Symbolics.array_term(convert, ptype, chain_params, size=size(chain_params))

lazyconvert_params = Symbolics.array_term((x, axes) -> ComponentVector(x, axes), chain_params, getaxes(init_params), size=size(chain_params))

@variables nn_in[1:2]
@variables nn_out[1:2]

expr = LuxCore.stateless_apply(chain, nn_in, lazyconvert_params)
ass = [Assignment(nn_in, MakeArray([a, b], Vector)), Assignment(nn_out, expr), Assignment(c, nn_out[1]), Assignment(d, nn_out[2])]
let_ = Let(ass, c + d, false)
func = @RuntimeGeneratedFunction(
    toexpr(Func([a, b, ptype, p], [], let_))
)
func(1, 2, typeof(ComponentVector(init_params)), collect(init_params))
