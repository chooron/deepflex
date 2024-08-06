# 导入模块
using ModelingToolkit
using CSV
using DataFrames
using ComponentArrays
using BenchmarkTools
using NamedTupleTools
using Plots
include("../../src/LumpedHydro.jl")

unit = LumpedHydro.ExpHydro.Unit(name=:exphydro)

f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = ComponentVector(snowwater=0.0, soilwater=1303.004248)
pas = ComponentVector(params=params, initstates=init_states)

file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:100)
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
solver = LumpedHydro.ODESolver()
result = unit(input, pas, timeidx=ts, solver=solver)

plot(result.soilwater)
# plot!(df[ts,"flow(mm)"])