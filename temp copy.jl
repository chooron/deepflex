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
    toexpr(Func([DestructuredArgs([a, b, c, d, e])], [], let_))
)

chain = Lux.Chain(Lux.Dense(2, 16), Lux.Dense(16, 2))
init_params = Lux.initialparameters(StableRNG(42), chain)
chain_params = first(@parameters p[1:length(init_params)] = Vector(ComponentVector(init_params)))
@parameters ptype::typeof(typeof(init_params)) = typeof(init_params) [tunable = false]
lazyconvert_params = Symbolics.array_term(convert, ptype, chain_params, size=size(chain_params))
@variables nn_in[1:2]
@variables nn_out[1:2]

expr = LuxCore.stateless_apply(chain, nn_in, lazyconvert_params)
ass = [Assignment(nn_out, expr)]
let_ = Let(ass, nn_out, false)
func = @RuntimeGeneratedFunction(
    toexpr(Func([nn_in, ptype, p], [], let_))
)
nnf = (x, p) -> LuxCore.stateless_apply(chain, x, p)
nnf([1, 2], init_params)

func([1, 2], typeof((init_params)), (init_params))
