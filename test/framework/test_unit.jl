# 导入模块
using ModelingToolkit
using CSV
using DataFrames
using ComponentArrays
using BenchmarkTools
using NamedTupleTools
using DataFrames
include("../../src/DeepFlex.jl")

model = DeepFlex.ExpHydro(name=:sf)

f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
params = (f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = (snowwater=0.0, soilwater=1303.004248)

file_path = "data/camels/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = 1:10000
input = (time=ts, lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
@benchmark results = model(input, params, init_states) # 分开计算是55.3ms
# @btime sol = DeepFlex.solve_prob(ele, input=input, params=params, init_states=init_states) # 1.47s
# @btime sol = DeepFlex.solve_probv2(ele, input=input, params=params, init_states=init_states) # 35.292ms