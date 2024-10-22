# 导入模块
using CSV
using DataFrames
using ComponentArrays
using BenchmarkTools
using NamedTupleTools

include("../../../src/HydroModels.jl")

ele = HydroModels.ExpHydro.SurfaceStorage(name=:sf)

f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
ps = [f, Smax, Qmax, Df, Tmax, Tmin]
params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = ComponentVector(snowpack=0.0, soilwater=1303.004248)
pas = ComponentVector(params=params, initstates=init_states)

file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
input_matrix = Matrix(reduce(hcat, [input[1], input[2], input[3]])') # (var nm * ts len)
solver = HydroModels.ODESolver()
config = (solver=solver, )
results = ele(input_matrix, pas, ts, config=config, convert_to_ntp=false)
