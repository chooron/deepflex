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
include("../models/m50.jl")

# load data
data = CSV.read("data/exphydro/01013500.csv", DataFrame)

# predefine the parameters
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

# load data
file_path = "data/m50/01013500.csv"
data = CSV.File(file_path)
df = DataFrame(data)
ts = collect(1:1000)
# cols: Baseflow,Evap,Flow,Infiltration,Lday,Melt,Pet,Prcp,Rainfall,Snowfall,Surfaceflow,Temp,SoilWater,SnowWater
lday_vec = df[ts, "Lday"]
prcp_vec = df[ts, "Prcp"]
temp_vec = df[ts, "Temp"]
flow_vec = df[ts, "Flow"]

log_flow_vec = log.(flow_vec)
log_evap_div_lday_vec = log.(df[ts, "Evap"] ./ lday_vec)
norm_prcp_vec = (prcp_vec .- mean(prcp_vec)) ./ std(prcp_vec)
norm_temp_vec = (temp_vec .- mean(temp_vec)) ./ std(temp_vec)
norm_snw_vec = (df[ts, "SnowWater"] .- mean(df[ts, "SnowWater"])) ./ std(df[ts, "SnowWater"])
norm_slw_vec = (df[ts, "SoilWater"] .- mean(df[ts, "SoilWater"])) ./ std(df[ts, "SoilWater"])
nn_input = (norm_snw=norm_snw_vec, norm_slw=norm_slw_vec, norm_temp=norm_temp_vec, norm_prcp=norm_prcp_vec)

ep_grad_opt = HydroModels.GradOptimizer(component=ep_nn_flux, solve_alg=Adam(1e-2), adtype=Optimization.AutoZygote(), maxiters=1000)
ep_input_matrix = Matrix(reduce(hcat, collect(nn_input[HydroModels.get_input_names(ep_nn_flux)]))')
ep_output = (log_evap_div_lday=log_evap_div_lday_vec,)

ep_opt_params, epnn_loss_df = ep_grad_opt(
    [ep_input_matrix], [ep_output],
    tunable_pas=ComponentVector(nn=(epnn=ep_nn_params,)),
    const_pas=ComponentVector(),
    return_loss_df=true
)

q_grad_opt = HydroModels.GradOptimizer(component=q_nn_flux, solve_alg=Adam(1e-2), adtype=Optimization.AutoZygote(), maxiters=1000)
q_input_matrix = Matrix(reduce(hcat, collect(nn_input[HydroModels.get_input_names(q_nn_flux)]))')
q_output = (log_flow=log_flow_vec,)

q_opt_params, qnn_loss_df = q_grad_opt(
    [q_input_matrix], [q_output],
    tunable_pas=ComponentVector(nn=(qnn=q_nn_params,)),
    const_pas=ComponentVector(),
    return_loss_df=true
)

m50_opt = HydroModels.GradOptimizer(component=m50_model, solve_alg=Adam(1e-2), adtype=Optimization.AutoForwardDiff(), maxiters=100)
# config = (solver=HydroModels.ODESolver(sensealg=BacksolveAdjoint(autodiff=ZygoteVJP())), interp=LinearInterpolation)
# config = (solver=HydroModels.DiscreteSolver(sensealg=ZygoteAdjoint()),)
norm_pas = ComponentVector(
    snowpack_mean=mean(norm_snw_vec), soilwater_mean=mean(norm_slw_vec), prcp_mean=mean(norm_prcp_vec), temp_mean=mean(norm_temp_vec),
    snowpack_std=std(norm_snw_vec), soilwater_std=std(norm_slw_vec), prcp_std=std(norm_prcp_vec), temp_std=std(norm_temp_vec)
)
m50_const_pas = ComponentVector(
    initstates=ComponentVector(snowpack=0.0, soilwater=1300.0),
    params=ComponentVector(Tmin=Tmin, Tmax=Tmax, Df=Df; norm_pas...)
)
m50_tunable_pas = ComponentVector(
    nn=ComponentVector(epnn=ep_nn_params, qnn=q_nn_params)
)
m50_input = (prcp=prcp_vec, lday=lday_vec, temp=temp_vec)
pas = ComponentVector(
    initstates=ComponentVector(snowpack=0.0, soilwater=1300.0),
    params=ComponentVector(Tmin=Tmin, Tmax=Tmax, Df=Df; norm_pas...),
    nn=ComponentVector(epnn=ep_nn_params, qnn=q_nn_params)
)
HydroModels.update_ca(pas, m50_tunable_pas)

# m50_opt_params, m50_loss_df = m50_opt(
#     [m50_input], [q_output],
#     tunable_pas=m50_tunable_pas,
#     const_pas=m50_const_pas,
#     config=[config],
#     return_loss_df=true
# )