using DataFrames
using CSV
using JLD2
using ComponentArrays
using HydroModels
using Plots
using OrdinaryDiffEq
using NamedTupleTools
include("../models/HBV.jl")

# 加载流域连接关系信息
grid_ids = collect(1:length(index_info))
param_types = Symbol.(grid_ids)
grid_area = (0.1 * 111)^2

timeidx = collect(1:length(qobs))[1:9000]
# 准备模型的输入数据:Vector{NamedTuple}
input_data = Vector{NamedTuple}()
for id in grid_ids
    push!(input_data, (prcp=grid_prec_df[timeidx, string(id)] .* 24, ep=grid_prec_df[timeidx, string(id)] .* 2.4))
end

# lumped model下的最优参数:  102.10581100749681  -2.57143626082281  61.515655243008084   2.3368971035090835
# 初始化模型参数, 假设每个子流域享有一个参数
basic_params = (x1=102.10581100749681, x2=-2.57143626082281, x3=61.515655243008084, x4=2.3368971035090835, lag=0.1)
initstates = (soilwater=0.0, routingstore=0.0, s_river=0.0)
area_coefs = NamedTuple{Tuple(param_types)}([(area_coef=grid_area,) for i in eachindex(param_types)])
node_params = NamedTuple{Tuple(param_types)}(repeat([basic_params], length(grid_ids)))
node_params = merge_recursive(node_params, area_coefs)
node_initstates = NamedTuple{Tuple(param_types)}(repeat([initstates], length(grid_ids)))
node_pas = ComponentVector(params=node_params, initstates=node_initstates)

gr4j_grid_models = load_gr4j_grid_model(flwdir_matrix, index_info, param_types)
base_config = (solver=ODESolver(),)
timeidx_hourly = collect(timeidx .- 1) .* 12
# route_config = (timeidx=timeidx_hourly, solver=ODESolver(alg=Rosenbrock23(), reltol=1e-2, abstol=1e-2),)
route_config = (timeidx=timeidx_hourly, solver=ODESolver(),)
configs = [base_config, base_config, NamedTuple(), route_config] 

n = HydroModels.build_grid_digraph(flwdir_matrix, index_info)

# 运行模型
output = gr4j_grid_models(input_data, node_pas, config=configs, convert_to_ntp=true);
# output_105 = output[105]
# # 获取输出节点数据105
# outlet_result = DataFrame(output_105)
# # outlet_result = DataFrame(output[1])
# # node_2_result = DataFrame(output[2])
# # # 绘制结果
# plot(outlet_result.q_routed[1000:1500], label="Qsim_routed")
# plot!(qobs[1000:1500], label="Qobs")