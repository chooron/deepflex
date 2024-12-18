using Lux
using ComponentArrays
using Symbolics
using SymbolicUtils
using SymbolicUtils.Code
using StableRNGs
using ModelingToolkit
using BenchmarkTools

q_nn = Lux.Chain(
    Lux.Dense(2 => 16, tanh),
    Lux.Dense(16 => 16, leakyrelu),
    Lux.Dense(16 => 16, leakyrelu),
    Lux.Dense(16 => 1, leakyrelu),
    name=:qnn
)
q_nn_ps, q_nn_st = Lux.setup(StableRNGs.LehmerRNG(1234), q_nn)
q_NN_stateful = Lux.StatefulLuxLayer{true}(q_nn, nothing, q_nn_st) # 不能用symbolic表示
q_nn_params_ca = ComponentVector(first(Lux.setup(StableRNGs.LehmerRNG(1234), q_nn)))
q_nn_states_ca = ComponentVector((Lux.setup(StableRNGs.LehmerRNG(1234), q_nn))[2])
q_nn_params_vec = Vector(q_nn_params_ca)
q_nn_states_vec = Vector(q_nn_states_ca)
q_axes = getaxes(q_nn_params_ca)
q_nn_func1 = (x, q) -> LuxCore.stateless_apply(q_nn, x, ComponentVector(q, q_axes))
q_nn_func2 = (x, q) -> LuxCore.stateless_apply(q_nn, x, q)
q_nn_func3 = (x, q) -> LuxCore.apply(q_NN_stateful, x, q)

chain_params = first(@parameters qnn_ps[1:length(q_nn_params_ca)] = Vector(q_nn_params_ca))
lazy_params = Symbolics.array_term((x, axes) -> ComponentVector(x, axes), chain_params, q_axes, size=size(chain_params))

nn_input = first(@variables x[1:2])
nn_output = first(@variables y[1:1])
nn_input_vars = @variables a b c
flux_expr = LuxCore.apply(q_nn, nn_input, lazy_params, q_nn_st)
flux_expr2 = LuxCore.stateless_apply(q_nn, nn_input, lazy_params)

assign_list = [
    Assignment(nn_input, MakeArray([a, b], Vector)),
    Assignment(nn_output, flux_expr[1]),
    Assignment(c, nn_output[1]),
]
outputs_arr = MakeArray([c], Vector)
func_args = [DestructuredArgs([a, b]), DestructuredArgs([chain_params])]

call_func = eval(toexpr(Func(func_args, [], Let(assign_list, outputs_arr, false))))

input_matrix = rand(2, 1000)
q_nn_func3 = (x) -> q_nn_func2(x, ComponentVector(q_nn_params_ca, q_axes))
@btime q_nn_func1.(eachslice(input_matrix, dims=2), Ref(q_nn_params_vec))
@btime q_nn_func2.(eachslice(input_matrix, dims=2), Ref(q_nn_params_ca))
@btime q_nn_func3.(eachslice(input_matrix, dims=2))
# @btime call_func.(eachslice(input_matrix, dims=2), Ref([q_nn_params_vec]))
# @btime [ComponentVector(q_nn_params_vec, q_axes) for _ in 1:1000];