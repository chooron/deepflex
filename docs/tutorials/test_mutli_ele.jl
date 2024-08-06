#* 测试多个element的并行求解
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

file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)

node_names = [Symbol(:node, i) for i in 1:100]
ele = LumpedHydro.ExpHydro.SurfaceStorage(name=:sf)
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = ComponentVector(snowwater=0.0, soilwater=1303.004248)

pas = ComponentVector(params=params, initstates=init_states)
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])

node_pas = ComponentVector(NamedTuple{Tuple(node_names)}(repeat([pas], length(node_names))))
node_params = ComponentVector(NamedTuple{Tuple(node_names)}([node_pas[nm][:params] for nm in node_names]))
node_initstates = ComponentVector(NamedTuple{Tuple(node_names)}([node_pas[nm][:initstates] for nm in node_names]))
node_input = NamedTuple{Tuple(node_names)}(repeat([input], length(node_names)))

result = LumpedHydro.solve_multi_prob(ele, input=node_input, params=node_params, init_states=node_initstates, timeidx=ts)
