# import lib
using CSV
using DataFrames
using ComponentArrays
using StructArrays
using OptimizationOptimisers
using BenchmarkTools
using NamedTupleTools
using Lux
using StableRNGs
using Enzyme
using OrdinaryDiffEq
include("../../src/LumpedHydro.jl")

# test exphydro model

# base param names
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

# input data
file_path = "data/m50/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)
input = (m50=StructArray(lday=df[ts, "Lday"], temp=df[ts, "Temp"], prcp=df[ts, "Prcp"], snowwater=df[ts, "SnowWater"], infiltration=df[ts, "Infiltration"]),)
output = (flow=df[ts, "Flow"],)

mean_snowwater, std_snowwater = mean(df[ts, "SnowWater"]), std(df[ts, "SnowWater"])
mean_soilwater, std_soilwater = mean(df[ts, "SoilWater"]), std(df[ts, "SoilWater"])
mean_temp, std_temp = mean(df[ts, "Temp"]), std(df[ts, "Temp"])
mean_prcp, std_prcp = mean(df[ts, "Prcp"]), std(df[ts, "Prcp"])

et_ann = Lux.Chain(Lux.Dense(3, 16, Lux.tanh), Lux.Dense(16, 1, Lux.leakyrelu)) # Lux.Dense(16, 16, Lux.leakyrelu), 
q_ann = Lux.Chain(Lux.Dense(2, 16, Lux.tanh), Lux.Dense(16, 1, Lux.leakyrelu))  # Lux.Dense(16, 16, Lux.leakyrelu), 
et_ann_p = LuxCore.initialparameters(StableRNG(42), et_ann)
q_ann_p = LuxCore.initialparameters(StableRNG(42), q_ann)

tunable_pas = ComponentVector(m50=(params=ComponentVector(
    f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin,
    etnn=et_ann_p, qnn=q_ann_p),),)

const_pas = ComponentVector(m50=(params=ComponentVector(
        mean_snowwater=mean_snowwater, std_snowwater=std_snowwater,
        mean_soilwater=mean_soilwater, std_soilwater=std_soilwater,
        mean_temp=mean_temp, std_temp=std_temp,
        mean_prcp=mean_prcp, std_prcp=std_prcp), initstates=ComponentVector(snowwater=0.0, soilwater=1300.0), weight=1.0),)

params_axes = getaxes(tunable_pas)

model = LumpedHydro.M50.Node(name=:m50, mtk=false, step=false)

best_pas = LumpedHydro.param_grad_optim(
    model,
    tunable_pas=tunable_pas,
    const_pas=const_pas,
    input=input,
    target=output,
    timeidx=ts,
    adtype=Optimization.AutoDiffractor(),
    solver=LumpedHydro.ODESolver(alg=Rosenbrock23()),
    maxiters=100
)

# total_params = LumpedHydro.merge_ca(best_pas, const_pas)[:param]
# result = model(input, total_params, step=false)