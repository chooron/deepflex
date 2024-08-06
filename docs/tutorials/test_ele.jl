# 导入模块
using CSV
using DataFrames
using ComponentArrays

using BenchmarkTools
using NamedTupleTools
using OrdinaryDiffEq
using Zygote
using ModelingToolkit
using Lux

include("../../src/LumpedHydro.jl")

ele = LumpedHydro.ExpHydro.SurfaceStorage(name=:sf)

f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
ps = [f, Smax, Qmax, Df, Tmax, Tmin]
params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = ComponentVector(snowwater=0.0, soilwater=1303.004248)
pas = ComponentVector(params=params, initstates=init_states)

file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
input_matrix = reduce(hcat, [input[1],input[2],input[3]])' # (var nm * ts len)
solver = LumpedHydro.ODESolver()
results = ele(input_matrix, pas, timeidx=ts, solver=solver)
# @btime [LumpedHydro.solve_single_prob(ele, input=input, params=params, init_states=init_states, timeidx=ts) for _ in 1:100]