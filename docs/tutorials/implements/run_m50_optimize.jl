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
using JLD2
using OptimizationOptimisers
using SciMLSensitivity
using HydroModels
# include("../../../src/HydroModels.jl")
include("../models/m50.jl")

# load data
data = CSV.read("../data/exphydro/01013500.csv", DataFrame)
ts = collect(1:10000)
input = (lday=data[ts, "dayl(day)"], temp=data[ts, "tmean(C)"], prcp=data[ts, "prcp(mm/day)"])
flow_vec = data[ts, "flow(mm)"]
load_exphydro = load("tutorials/implements/save/exphydro_opt.jld2")
exphydro_output = load_exphydro["output"]
exphydro_opt_params = load_exphydro["opt_params"]
# predefine the parameters
Df, Tmax, Tmin = exphydro_opt_params.params.Df, exphydro_opt_params.params.Tmax, exphydro_opt_params.params.Tmin

#* prepare the nn input
log_flow_vec = log.(exphydro_output.flow)
log_evap_div_lday_vec = log.(exphydro_output.evap ./ exphydro_output.lday)
norm_prcp_vec = (exphydro_output.prcp .- mean(exphydro_output.prcp)) ./ std(exphydro_output.prcp)
norm_temp_vec = (exphydro_output.temp .- mean(exphydro_output.temp)) ./ std(exphydro_output.temp)
norm_snw_vec = (exphydro_output.snowpack .- mean(exphydro_output.snowpack)) ./ std(exphydro_output.snowpack)
norm_slw_vec = (exphydro_output.soilwater .- mean(exphydro_output.soilwater)) ./ std(exphydro_output.soilwater)
nn_input = (norm_snw=norm_snw_vec, norm_slw=norm_slw_vec, norm_temp=norm_temp_vec, norm_prcp=norm_prcp_vec)
ep_input_matrix = Matrix(reduce(hcat, collect(nn_input[HydroModels.get_input_names(ep_nn_flux)]))')
q_input_matrix = Matrix(reduce(hcat, collect(nn_input[HydroModels.get_input_names(q_nn_flux)]))')

#* train the ep nn
ep_grad_opt = HydroModels.GradOptimizer(component=ep_nn_flux, solve_alg=Adam(1e-2), adtype=AutoZygote(), maxiters=1000)
ep_output = (log_evap_div_lday=log_evap_div_lday_vec,)
ep_opt_params, epnn_loss_df = ep_grad_opt(
    [ep_input_matrix], [ep_output],
    tunable_pas=ComponentVector(nn=(epnn=ep_nn_params,)),
    const_pas=ComponentVector(),
    return_loss_df=true
)
ep_flux_output = ep_nn_flux(ep_input_matrix, ep_opt_params)
#* train the q nn
q_grad_opt = HydroModels.GradOptimizer(component=q_nn_flux, solve_alg=Adam(1e-2), adtype=AutoZygote(), maxiters=1000)
q_opt_params, qnn_loss_df = q_grad_opt(
    [q_input_matrix], [(log_flow=log_flow_vec,)],
    tunable_pas=ComponentVector(nn=(qnn=q_nn_params,)),
    const_pas=ComponentVector(),
    return_loss_df=true
)
q_flux_output = q_nn_flux(q_input_matrix, q_opt_params)

# save("doc/tutorials/implements/save/m50_nn_opt.jld2", "epnn_params", ep_opt_params, "qnn_params", q_opt_params, "epnn_loss_df", epnn_loss_df, "qnn_loss_df", qnn_loss_df)

#* define parameters
norm_pas = ComponentVector(
    snowpack_mean=mean(exphydro_output.snowpack), soilwater_mean=mean(exphydro_output.soilwater), prcp_mean=mean(input.prcp), temp_mean=mean(input.temp),
    snowpack_std=std(exphydro_output.snowpack), soilwater_std=std(exphydro_output.soilwater), prcp_std=std(input.prcp), temp_std=std(input.temp)
)
m50_const_pas = ComponentVector(
    initstates=ComponentVector(snowpack=0.0, soilwater=1300.0),
    params=ComponentVector(Tmin=Tmin, Tmax=Tmax, Df=Df; norm_pas...)
)
m50_tunable_pas = ComponentVector(
    nn=ComponentVector(epnn=ep_nn_params, qnn=q_nn_params)
)
#* build optimizer for m50 model
m50_opt = HydroModels.GradOptimizer(component=m50_model, solve_alg=Adam(1e-2), adtype=AutoZygote(), maxiters=100)

#* prepare input and config
m50_input = (prcp=exphydro_output.prcp, lday=exphydro_output.lday, temp=exphydro_output.temp)
config = (solver=HydroModels.ODESolver(sensealg=BacksolveAdjoint(autojacvec=ZygoteVJP())), interp=LinearInterpolation, alg=BS3())
q_output = (flow=flow_vec,)

#* optimize the model
m50_opt_params, m50_loss_df = m50_opt(
    [m50_input], [q_output],
    tunable_pas=m50_tunable_pas,
    const_pas=m50_const_pas,
    config=[config],
    return_loss_df=true
)
save("doc/tutorials/implements/save/m50_opt.jld2", "loss_df", m50_loss_df, "opt_params",
    m50_opt_params, "epnn_params", ep_opt_params, "qnn_params", q_opt_params)
plot(m50_loss_df.loss)
m50_opt_params = load("doc/tutorials/implements/save/m50_opt.jld2")["opt_params"]
m50_output = m50_model(input, m50_opt_params, config=config, convert_to_ntp=true)

tmp_params = ComponentVector(
    initstates=ComponentVector(snowpack=0.0, soilwater=1300.0),
    params=ComponentVector(Tmin=Tmin, Tmax=Tmax, Df=Df; norm_pas...),
    nn=ComponentVector(epnn=ep_opt_params, qnn=q_opt_params)
)
m50_outputv2 = m50_model(input, tmp_params, config=config, convert_to_ntp=true)

result = (exphydro=exphydro_output.flow, m50=m50_output.flow, qnn=exp.((q_flux_output)[1, :]), sim=flow_vec)
save("doc/tutorials/implements/save/m50_result.jld2",
    "exphydro", exphydro_output.flow, "m50", m50_output.flow,
    "qnn", exp.((q_flux_output)[1, :]), "obs", flow_vec)