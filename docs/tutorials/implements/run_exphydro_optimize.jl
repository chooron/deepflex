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
using OptimizationBBO
using SciMLSensitivity
using JLD2
# using HydroModels
include("../../../src/HydroModels.jl")
include("../models/exphydro.jl")
# model_grad_opt = HydroModels.GradOptimizer(component=exphydro_model, solve_alg=Adam(1e-2), adtype=Optimization.AutoForwardDiff(), maxiters=100)


# load data
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
q_vec = df[ts, "flow(mm)"]

model_hydro_opt = HydroModels.HydroOptimizer(
    component=exphydro_model,
    maxiters=10000,
    warmup=100,
    solve_alg=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    loss_func=(obs, sim) -> sum((obs .- sim) .^ 2) / length(obs)
)
tunable_pas = ComponentVector(params=ComponentVector(
    f=0.01674478, Smax=1709.461015, Qmax=18.46996175,
    Df=2.674548848, Tmax=0.175739196, Tmin=-2.092959084
))
const_pas = ComponentVector(initstates=ComponentVector(snowpack=0.0, soilwater=1300.0))
config = (solver=HydroModels.ODESolver(), interp=LinearInterpolation)

exphydro_hydro_opt_params, loss_df = model_hydro_opt(
    [input], [(flow=q_vec,)],
    tunable_pas=tunable_pas,
    const_pas=const_pas,
    config=[config],
    lb=[0.0, 100.0, 10.0, 0.0, 0.0, -3.0],
    ub=[0.1, 2000.0, 50.0, 5.0, 3.0, 0.0],
    return_loss_df=true
)
output = exphydro_model(input, exphydro_hydro_opt_params, config=config, convert_to_ntp=true)
save("doc/tutorials/implements/save/exphydro_opt.jld2", "loss_df", loss_df, "opt_params", exphydro_hydro_opt_params, "output", output)
# exphydro_grad_opt_params, loss_df = model_grad_opt(
#     [input], [(flow=q_vec,)],
#     tunable_pas=tunable_pas,
#     const_pas=const_pas,
#     config=[config],
#     return_loss_df=true
# )