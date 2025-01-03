# 导入模块
using CSV
using DataFrames
using ComponentArrays
using ModelingToolkit
using Lux
using BenchmarkTools
using DataInterpolations
using CUDA
using Zygote

include("../src/HydroModels.jl")

HydroFlux = HydroModels.HydroFlux
StateFlux = HydroModels.StateFlux
HydroBucket = HydroModels.HydroBucket
HydroModel = HydroModels.HydroModel
include("../models/exphydro.jl")

ele = bucket_1
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
# 0.0167, 1709.46, 18.47, 2.6745, 0.1757, -2.093
params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = ComponentVector(snowpack=0.0, soilwater=1303.004248)
pas = ComponentVector(params=params, initstates=init_states)
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:100)

# single node input
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"]) 
input_arr = Matrix(reduce(hcat, collect(input[ele.meta.inputs]))')
config = (solver=HydroModels.ManualSolver{false}(), timeidx=ts)
results = ele(input_arr, pas)


gradient((p) -> sum(ele(input_arr, p, config=config)), pas)
