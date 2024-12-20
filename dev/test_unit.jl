# 导入模块
using CSV
using DataFrames
using ComponentArrays
using HydroModelTools
using ModelingToolkit
using Plots

include("../../src/HydroModels.jl")
include("../src/models/exphydro.jl")

# define parameters and initial states
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = ComponentVector(snowpack=0.0, soilwater=1303.004248)
pas = ComponentVector(params=params, initstates=init_states)

# load data
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])

# run model with single node input
solver = ODESolver(reltol=1e-3, abstol=1e-3)
adapter = HydroModels.NamedTupleIOAdapter(exphydro_model)
result = adapter(input, pas, config=(solver=solver, timeidx=ts))

plot(result.flow)
plot!(df[!, "flow(mm)"])
# # # run model with multi node input
# node_num = 10
# inputs = repeat([input], node_num)
# ptypes = [Symbol(:node, i) for i in 1:node_num]
# params_multi = ComponentVector(NamedTuple{Tuple(ptypes)}(repeat([params], node_num)))
# init_states_multi = ComponentVector(NamedTuple{Tuple(ptypes)}(repeat([init_states], node_num)))
# pas_multi = ComponentVector(params=params_multi, initstates=init_states_multi)
# results = exphydro_model(inputs, pas_multi, config=(solver=solver, timeidx=ts), convert_to_ntp=true)
