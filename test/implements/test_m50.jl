# 导入模块
using CSV
using DataFrames
using ComponentArrays
using BenchmarkTools
using NamedTupleTools
using DataFrames
using Lux, LuxCore
using StableRNGs
using OrdinaryDiffEq
using ModelingToolkit
# using LumpedHydro

include("../../src/LumpedHydro.jl")
model = LumpedHydro.M50.Node(name=:m50, mtk=false, step=false)
# unit_sys = model.units[1]
# unknowns(unit_sys.system)
# base param names
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

# input data
file_path = "data/m50/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = 1:1000

mean_snowwater, std_snowwater = mean(df[ts, "SnowWater"]), std(df[ts, "SnowWater"])
mean_soilwater, std_soilwater = mean(df[ts, "SoilWater"]), std(df[ts, "SoilWater"])
mean_temp, std_temp = mean(df[ts, "Temp"]), std(df[ts, "Temp"])
mean_prcp, std_prcp = mean(df[ts, "Prcp"]), std(df[ts, "Prcp"])

input = (m50=(time=ts, lday=df[ts, "Lday"], temp=df[ts, "Temp"], prcp=df[ts, "Prcp"]),) # , snowwater=df[ts, "SnowWater"], infiltration=df[ts, "Infiltration"]),

et_ann = Lux.Chain(Lux.Dense(3, 16, Lux.tanh), Lux.Dense(16, 1, Lux.leakyrelu))
q_ann = Lux.Chain(Lux.Dense(2, 16, Lux.tanh), Lux.Dense(16, 1, Lux.leakyrelu))
et_ann_p = LuxCore.initialparameters(StableRNG(42), et_ann)
q_ann_p = LuxCore.initialparameters(StableRNG(42), q_ann)

params = ComponentVector(
    f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin,
    mean_snowwater=mean_snowwater, std_snowwater=std_snowwater,
    mean_soilwater=mean_soilwater, std_soilwater=std_soilwater,
    mean_temp=mean_temp, std_temp=std_temp,
    mean_prcp=mean_prcp, std_prcp=std_prcp,
    etnn=et_ann_p, qnn=q_ann_p)
initstates = ComponentVector(snowwater=0.0, soilwater=1303.004248)
pas = ComponentVector(m50=(params=params, initstates=initstates, weight=1.0),)

solver = LumpedHydro.ODESolver(alg=Rosenbrock23())
@btime results = model(input, pas, solver=solver)

# q_ann = Lux.Chain(
#     Lux.Dense(2 => 16, Lux.tanh),
#     # Lux.Dense(16 => 16, Lux.leakyrelu),
#     Lux.Dense(16 => 1, Lux.leakyrelu)
# )
# func = (x, p) -> LuxCore.stateless_apply(q_ann, x, p)

