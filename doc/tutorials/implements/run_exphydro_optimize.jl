using CSV
using DataFrames
using Lux
using ModelingToolkit
using LuxCore
using StableRNGs
using ComponentArrays
using DataInterpolations
using OrdinaryDiffEq
using Statistics
using BenchmarkTools
using Plots
using OptimizationOptimisers
using SciMLSensitivity
# using HydroModels
include("../../../src/HydroModels.jl")
include("../models/exphydro.jl")

# load data
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
q_vec = df[ts, "flow(mm)"]

model_grad_opt = HydroModels.GradOptimizer(component=model, solve_alg=Adam(1e-2), adtype=Optimization.AutoForwardDiff(), maxiters=100)
model_hydro_opt = HydroModels.HydroOptimizer(component=model, maxiters=100)
tunable_params = ComponentVector(f=0.01674478, Smax=1709.461015, Qmax=18.46996175, Df=2.674548848, Tmax=0.175739196, Tmin=-2.092959084)
tunable_pas = ComponentVector(params=tunable_params)
const_pas = ComponentVector(initstates=ComponentVector(snowpack=0.0, soilwater=1300.0))
config = (solver=HydroModels.ODESolver(), interp=LinearInterpolation)

exphydro_hydro_opt_params, loss_dfv2 = model_hydro_opt(
    [input], [(flow=q_vec,)],
    tunable_pas=tunable_pas,
    const_pas=const_pas,
    config=[config],
    lb=[0.0, 100.0, 10.0, 0.0, 0.0, -3.0],
    ub=[0.1, 2000.0, 50.0, 5.0, 3.0, 0.0],
    return_loss_df=true
)

exphydro_grad_opt_params, loss_df = model_grad_opt(
    [input], [(flow=q_vec,)],
    tunable_pas=tunable_pas,
    const_pas=const_pas,
    config=[config],
    return_loss_df=true
)