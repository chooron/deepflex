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


# include need add the module name: HydroModels
HydroFlux = HydroModels.HydroFlux
StateFlux = HydroModels.StateFlux
NeuralFlux = HydroModels.NeuralFlux
HydroBucket = HydroModels.HydroBucket
HydroModel = HydroModels.HydroModel
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

#! parameters in the Exp-Hydro model
@parameters Tmin Tmax Df Smax f Qmax
#! parameters in normalize flux
@parameters snowpack_std snowpack_mean
@parameters soilwater_std soilwater_mean
@parameters prcp_std prcp_mean
@parameters temp_std temp_mean

#! hydrological flux in the Exp-Hydro model
@variables prcp temp lday pet rainfall snowfall
@variables snowpack soilwater lday pet
@variables melt log_evap_div_lday log_flow asinh_melt asinh_ps asinh_pr
@variables norm_snw norm_slw norm_temp norm_prcp

#! define the m100 NN
m100_nn = Lux.Chain(
    Lux.Dense(4, 32, tanh),
    Lux.Dense(32, 32, leakyrelu),
    Lux.Dense(32, 32, leakyrelu),
    Lux.Dense(32, 32, leakyrelu),
    Lux.Dense(32, 32, leakyrelu),
    Lux.Dense(32, 5),
    name=:m100nn
)
m100_nn_params = Vector(ComponentVector(first(Lux.setup(StableRNGs.LehmerRNG(1234), m100_nn))))

#! get init parameters for each NN

#! define the soil water reservoir
m100_funcs = [
    #* normalize
    HydroFlux([snowpack] => [norm_snw], [snowpack_mean, snowpack_std], exprs=[(snowpack - snowpack_mean) / snowpack_std]),
    HydroFlux([soilwater] => [norm_slw], [soilwater_mean, soilwater_std], exprs=[(soilwater - soilwater_mean) / soilwater_std]),
    HydroFlux([prcp] => [norm_prcp], [prcp_mean, prcp_std], exprs=[(prcp - prcp_mean) / prcp_std]),
    HydroFlux([temp] => [norm_temp], [temp_mean, temp_std], exprs=[(temp - temp_mean) / temp_std]),
    NeuralFlux([norm_snw, norm_slw, norm_prcp, norm_temp] => [log_evap_div_lday, log_flow, asinh_melt, asinh_ps, asinh_pr], m100_nn),
    HydroFlux([asinh_melt, snowpack] => [melt], exprs=[relu(sinh(asinh_melt) * step_func(snowpack))]),
]

m100_nn_flux = soil_funcs[end-1]
state_expr1 = relu(sinh(asinh_ps)) * step_func(-temp) - melt
state_expr2 = relu(sinh(asinh_pr)) + melt - step_func(soilwater) * lday * exp(log_evap_div_lday) - step_func(soilwater) * exp(log_flow)
m100_dfuncs = [
    StateFlux([asinh_ps, temp, melt], snowpack, Num[], expr=state_expr1),
    StateFlux([asinh_pr, melt, soilwater, lday, log_evap_div_lday, log_flow], soilwater, Num[], expr=state_expr2),
]
m100_bucket = HydroBucket(name=:m100_bucket, funcs=m100_funcs, dfuncs=m100_dfuncs)
m100_bucket.ode_func
#! define the Exp-Hydro model
m100_model = HydroModel(name=:m100, components=[m100_bucket]);

# load data
data = CSV.read("data/exphydro/01013500.csv", DataFrame)

# predefine the parameters
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

# load data
file_path = "data/m50/01013500.csv"
data = CSV.File(file_path)
df = DataFrame(data)
ts = collect(1:10000)
# cols: Baseflow,Evap,Flow,Infiltration,Lday,Melt,Pet,Prcp,Rainfall,Snowfall,Surfaceflow,Temp,SoilWater,SnowWater
lday_vec = df[ts, "Lday"]
prcp_vec = df[ts, "Prcp"]
temp_vec = df[ts, "Temp"]
flow_vec = df[ts, "Flow"]

log_flow_vec = log.(flow_vec)
log_evap_div_lday_vec = log.(df[ts, "Evap"] ./ lday_vec)
asinh_melt_vec = asinh.(df[ts, "Melt"])
asinh_ps_vec = asinh.(df[ts, "Snowfall"])
asinh_pr_vec = asinh.(df[ts, "Rainfall"])
norm_prcp_vec = (prcp_vec .- mean(prcp_vec)) ./ std(prcp_vec)
norm_temp_vec = (temp_vec .- mean(temp_vec)) ./ std(temp_vec)
norm_snw_vec = (df[ts, "SnowWater"] .- mean(df[ts, "SnowWater"])) ./ std(df[ts, "SnowWater"])
norm_slw_vec = (df[ts, "SoilWater"] .- mean(df[ts, "SoilWater"])) ./ std(df[ts, "SoilWater"])
nn_input = (norm_snw=norm_snw_vec, norm_slw=norm_slw_vec, norm_temp=norm_temp_vec, norm_prcp=norm_prcp_vec)

nn_grad_opt = HydroModels.GradOptimizer(component=m100_nn_flux, solve_alg=Adam(1e-2), adtype=Optimization.AutoZygote(), maxiters=1000)
nn_input_matrix = Matrix(reduce(hcat, collect(nn_input[HydroModels.get_input_names(m100_nn_flux)]))')
nn_output = (
    log_evap_div_lday=log_evap_div_lday_vec, log_flow=log_flow_vec,
    asinh_melt=asinh_melt_vec, asinh_ps=asinh_ps_vec, asinh_pr=asinh_pr_vec
)

nn_opt_params, nn_loss_df = nn_grad_opt(
    [nn_input_matrix], [nn_output],
    tunable_pas=ComponentVector(nn=(m100nn=m100_nn_params,)),
    const_pas=ComponentVector(),
    return_loss_df=true
)
norm_pas = ComponentVector(
    snowpack_mean=mean(norm_snw_vec), soilwater_mean=mean(norm_slw_vec), prcp_mean=mean(norm_prcp_vec), temp_mean=mean(norm_temp_vec),
    snowpack_std=std(norm_snw_vec), soilwater_std=std(norm_slw_vec), prcp_std=std(norm_prcp_vec), temp_std=std(norm_temp_vec)
)
m100_const_pas = ComponentVector(
    initstates=ComponentVector(snowpack=0.0, soilwater=1300.0),
    params=norm_pas
)
m100_input = (prcp=prcp_vec, lday=lday_vec, temp=temp_vec)
# 5.080 s (15939665 allocations: 16.79 GiB)
# @btime m100_model(m100_input, ComponentVector(nn=(m100nn=nn_opt_params,), params=norm_pas, initstates=ComponentVector(snowpack=0.0, soilwater=1300.0)))

m100_opt = HydroModels.GradOptimizer(component=m100_model, solve_alg=Adam(1e-2), adtype=Optimization.AutoZygote(), maxiters=100)
config = (solver=HydroModels.ODESolver(sensealg=BacksolveAdjoint(autodiff=ZygoteVJP())), interp=LinearInterpolation)
norm_pas = ComponentVector(
    snowpack_mean=mean(norm_snw_vec), soilwater_mean=mean(norm_slw_vec), prcp_mean=mean(norm_prcp_vec), temp_mean=mean(norm_temp_vec),
    snowpack_std=std(norm_snw_vec), soilwater_std=std(norm_slw_vec), prcp_std=std(norm_prcp_vec), temp_std=std(norm_temp_vec)
)
m100_const_pas = ComponentVector(
    initstates=ComponentVector(snowpack=0.0, soilwater=1300.0),
    params=norm_pas
)
m100_tunable_pas = ComponentVector(
    nn=ComponentVector(m100nn=m100_nn_params)
)
m100_input = (prcp=prcp_vec, lday=lday_vec, temp=temp_vec)
m100_opt_params, m100_loss_df = m100_opt(
    [m100_input], [nn_output],
    tunable_pas=m100_tunable_pas,
    const_pas=m100_const_pas,
    config=[config],
    return_loss_df=true
)
# # plot(m50_loss_df[!,:loss])
# @btime result = m50_model(m50_input, m50_opt_params; config)
