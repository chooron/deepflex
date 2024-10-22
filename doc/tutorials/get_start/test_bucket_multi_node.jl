#* 测试多个element的并行求解
using CSV
using DataFrames
using ComponentArrays
using BenchmarkTools
using NamedTupleTools
using DataInterpolations

include("../../../src/HydroModels.jl")

file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)

node_names = [Symbol(:node, i) for i in 1:100]
ele = HydroModels.ExpHydro.SurfaceStorage(name=:sf)
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = ComponentVector(snowpack=0.0, soilwater=1303.004248)

input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])

node_params = ComponentVector(NamedTuple{Tuple(node_names)}(repeat([params], length(node_names))))
node_initstates = ComponentVector(NamedTuple{Tuple(node_names)}(repeat([init_states], length(node_names))))
node_pas = ComponentVector(params=node_params, initstates=node_initstates)

input_arr = reduce(hcat, collect(input[HydroModels.get_input_names(ele)]))
node_input = reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([input_arr], length(node_names)))
node_input = permutedims(node_input, (2, 3, 1))
# result = HydroModels.solve_multi_prob(ele, input=node_input, pas=node_pas, timeidx=ts)
run_kwgs = (ptypes=node_names, interpolator=LinearInterpolation)

result = ele(node_input, node_pas, ts, kwargs=run_kwgs)
node_input = cat(node_input, result, dims=1)