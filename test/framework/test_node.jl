# 导入模块
using ModelingToolkit
using CSV
using DataFrames
using ComponentArrays
using BenchmarkTools
using NamedTupleTools
using DataFrames
include("../../src/DeepFlex.jl")

model = DeepFlex.ExpHydro.Node(name=:exphydro, mtk=true, step=false)

f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

file_path = "data/camels/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = 1:10000

unit_params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
unit_init_states = ComponentVector(snowwater=0.0, soilwater=1303.004248)
unit_input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"], time=ts)

# prob = DeepFlex.setup_input(model, input=unit_input, time=ts)
# new_prob = DeepFlex.setup_prob(model, prob, params=unit_params, init_states=unit_init_states)
unit_init_states = ComponentVector(snowwater=0.0, soilwater=1303.004248)

input = (exphydro=unit_input,)
pas = ComponentVector(exphydro=(params=unit_params, initstates=unit_init_states, weight=1.0))
results = model(input, pas)