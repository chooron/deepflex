using DataFrames
using CSV
using JLD2
using ComponentArrays
using Plots
using HydroModels
using OrdinaryDiffEq
using NamedTupleTools

include("load_data.jl")
include("models/gr4j.jl")

# 加载流域连接关系信息
dg8, hybas_id = load("data/jld2/sub_basin_network_8.jld2", "dg8", "hybas_id")

# 准备模型的输入数据:Vector{NamedTuple}
input_data = Vector{NamedTuple}()
timeidx = collect(1:length(qobs))
for id in hybas_id
    push!(input_data, (prcp=vector_prec_df[timeidx, string(id)] .* 24, ep=vector_prec_df[timeidx, string(id)] .* 2.4))
end

# lumped model下的最优参数:  102.10581100749681  -2.57143626082281  61.515655243008084   2.3368971035090835
# 初始化模型参数, 假设每个子流域享有一个参数 subasin_areas
basic_params = (x1=202.10581100749681, x2=-2.57143626082281, x3=31.515655243008084, x4=3.3368971035090835, lag=0.1)
initstates = (soilwater=0.0, routingstore=0.0, s_river=0.0)
param_types = Symbol.(hybas_id)
area_coefs = @. 24 * 3600 / (subasin_areas * 1e6) * 1e3
area_coefs = NamedTuple{Tuple(param_types)}([(area_coef=area_coefs[i],) for i in eachindex(area_coefs)])
node_params = NamedTuple{Tuple(param_types)}(repeat([basic_params], length(hybas_id)))
node_params = merge_recursive(node_params, area_coefs)
node_initstates = NamedTuple{Tuple(param_types)}(repeat([initstates], length(hybas_id)))
node_pas = ComponentVector(params=node_params, initstates=node_initstates)
gr4j_vector_models = load_gr4j_vector_model(dg8, param_types)
# 运行模型
timeidx_hourly = collect(timeidx) .* 12
base_config = (solver=ODESolver(),)
convert_config = NamedTuple()
route_config = (solver=ODESolver(), delta_t=12.0, timeidx=timeidx_hourly)
configs = [base_config, base_config, convert_config, route_config]
@time output = gr4j_vector_models(input_data, node_pas, config=configs, convert_to_ntp=true)
outlet_result = DataFrame(output[1])
node_2_result = DataFrame(output[2])

# # 绘制结果
plot(output[1].q_routed[2500:2700], label="Qsim")
plot!(qobs[2500:2700], label="Qobs")
# using HydroErrors
# HydroErrors.r2(outlet_result.flow_routed, qobs[timeidx])
