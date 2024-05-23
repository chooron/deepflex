# 导入模块
using ModelingToolkit
using CSV
using DataFrames
using ComponentArrays
using BenchmarkTools
using NamedTupleTools
using DataFrames
using StructArrays
using CairoMakie

include("../../src/LumpedHydro.jl")

model = LumpedHydro.ExpHydro.Node(name=:exphydro, mtk=true, step=false)

f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = 1:10000

unit_params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
unit_init_states = ComponentVector(snowwater=0.0, soilwater=1303.004248)
unit_input = StructArray(lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"], time=ts)

unit_init_states = ComponentVector(snowwater=0.0, soilwater=1303.004248)

input = (exphydro=unit_input, time=collect(ts))
pas = ComponentVector(exphydro=(params=unit_params, initstates=unit_init_states, weight=1.0))
results = model(input, pas)

fig = Figure(size=(400, 300))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
lines!(ax, ts, flow_vec, color=:red)
lines!(ax, ts, results[:flow], color=:blue)
fig