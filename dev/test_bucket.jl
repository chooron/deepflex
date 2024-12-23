# 导入模块
using CSV
using DataFrames
using ComponentArrays
using ModelingToolkit
using Lux
using BenchmarkTools
using HydroModelTools
using DataInterpolations
using Zygote
using SciMLSensitivity

include("../src/HydroModels.jl")

HydroFlux = HydroModels.HydroFlux
StateFlux = HydroModels.StateFlux
HydroBucket = HydroModels.HydroBucket
HydroModel = HydroModels.HydroModel
include("../models/exphydro.jl")

ele = bucket_1
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = ComponentVector(snowpack=0.0, soilwater=1303.004248)
pas = ComponentVector(params=params, initstates=init_states)

file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)

# single node input
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"]) 
# solver = HydroModels.ManualSolver{true}()
solver = ODESolver(sensealg=BacksolveAdjoint(autojacvec=ZygoteVJP()))
config = (solver=solver,)
input_arr = Matrix(reduce(hcat, collect(input[ele.meta.inputs]))')
results = ele(input_arr, pas, config=config)

Zygote.gradient((p) -> sum(ele(input_arr, p, config=config)[end, :]), pas)

# multi node input
# node_num = 10
# node_names = [Symbol(:node, i) for i in 1:node_num]
# node_params = ComponentVector(NamedTuple{Tuple(node_names)}(repeat([params], length(node_names))))
# node_initstates = ComponentVector(NamedTuple{Tuple(node_names)}(repeat([init_states], length(node_names))))
# node_pas = ComponentVector(params=node_params, initstates=node_initstates)

# input_arr = reduce(hcat, collect(input[HydroModels.get_input_names(ele)]))
# node_input = reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([input_arr], length(node_names)))
# node_input = permutedims(node_input, (2, 3, 1))
# run_kwgs = (ptypes=node_names, timeidx=ts)

# result = ele(node_input, node_pas, config=config, kwargs=run_kwgs)

# test_arr = ones(3, 3, 100)
# # node_input = cat(node_input, result, dims=1)