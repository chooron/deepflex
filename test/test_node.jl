# import lib
using CSV
using DataFrames
using CairoMakie
using ComponentArrays
using CairoMakie: Axis

# test exphydro model
include("../src/DeepFlex.jl")

# load data
file_path = "data/camels/01013500.csv"

data = CSV.File(file_path);
df = DataFrame(data);
lday_vec = df[1:100, "dayl(day)"]
prcp_vec = df[1:100, "prcp(mm/day)"]
temp_vec = df[1:100, "tmean(C)"]
flow_vec = df[1:100, "flow(mm)"]
snowwater_vec = zeros(size(lday_vec))
soilwater_vec = ones(size(lday_vec)) .* 1300

# build model
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
parameters = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = ComponentVector(SnowWater=0.0, SoilWater=1303.004248)
model = DeepFlex.M50(id="M50", parameters=parameters, init_states=init_states)
# todo 内部模型需要提供预训练功能
DeepFlex.pretrain(model, input=input, config=config)

# BacksolveAdjoint(autojacvec=ZygoteVJP())
inputs = ComponentVector(Prcp=prcp_vec, Lday=lday_vec, Temp=temp_vec)
# ! 当前求解结果只输出一个，考虑可能是模型没有训练过
result = DeepFlex.get_output(model, input=inputs, step=false, sensealg=DeepFlex.default_node_sensealg, reltol=1e-1, abstol=1e-1)

result_df = DataFrame(Dict(k => result[k] for k in keys(result)))
# plot result
fig = Figure(size=(400, 300))
ax = Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
x = range(1, 1000, length=1000)
lines!(ax, x, flow_vec, color=:red)
lines!(ax, x, result[:Flow], color=:blue)
fig

