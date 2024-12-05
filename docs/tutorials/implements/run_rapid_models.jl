using DataFrames
using CSV
using JLD2
using ComponentArrays
using Plots
using HydroModels
using OrdinaryDiffEq

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
# 初始化模型参数, 假设每个子流域享有一个参数
basic_params = (x1=102.10581100749681, x2=-2.57143626082281, x3=61.515655243008084, x4=2.3368971035090835, k=1.0, x=0.5)
initstates = (soilwater=0.0, routingstore=0.0)
param_types = Symbol.(hybas_id)
node_params = NamedTuple{Tuple(param_types)}(repeat([basic_params], length(hybas_id)))
node_initstates = NamedTuple{Tuple(param_types)}(repeat([initstates], length(hybas_id)))
node_pas = ComponentVector(params=node_params, initstates=node_initstates)
gr4j_vector_models = load_gr4j_rapid_model(dg8, subasin_areas)
# 运行模型
timeidx_hourly = collect(timeidx) .* 12
base_config = (solver=ODESolver(),)
route_config = (solver=DiscreteSolver(), delta_t=12.0, timeidx=timeidx_hourly)
configs = [base_config, base_config, route_config]
output = gr4j_vector_models(input_data, node_pas, config=configs, convert_to_ntp=true)
outlet_result = DataFrame(output[1])
node_2_result = DataFrame(output[2])

# # 绘制结果
plot(outlet_result.flow_routed[9100:10000], label="Qsim")
plot!(qobs[9100:10000], label="Qobs")
# using HydroErrors
# HydroErrors.r2(outlet_result.flow_routed, qobs[timeidx])
