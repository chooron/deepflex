using CSV
using Lux
using LuxCore
using Random
using DataFrames
using Symbolics
using ComponentArrays
using OrdinaryDiffEq
using ModelingToolkit
using BenchmarkTools
using StableRNGs
using Optimization
using NamedTupleTools
using HydroErrors
include("../../src/HydroModels.jl")

# load data
df = DataFrame(CSV.File("data/m50/01013500.csv"));
ts = collect(1:10000)
prcp_vec = df[ts, "Prcp"]
temp_vec = df[ts, "Temp"]
dayl_vec = df[ts, "Lday"]
snowpack_vec = df[ts, "SnowWater"]
soilwater_vec = df[ts, "SoilWater"]
qobs_vec = df[ts, "Flow"]

inputs = [prcp_vec, temp_vec, snowpack_vec, soilwater_vec]
means = mean.(inputs)
stds = std.(inputs)
(prcp_norm_vec, temp_norm_vec, snowpack_norm_vec, soilwater_norm_vec) = [
    @.((tmp_vec - mean) / std) for (tmp_vec, mean, std) in zip(inputs, means, stds)
]

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
@variables melt log_evap_div_lday log_flow
@variables norm_snw norm_slw norm_temp norm_prcp

SimpleFlux = HydroModels.SimpleFlux
StdMeanNormFlux = HydroModels.StdMeanNormFlux
NeuralFlux = HydroModels.NeuralFlux
StateFlux = HydroModels.StateFlux
HydroBucket = HydroModels.HydroBucket
HydroModel = HydroModels.HydroModel
step_func = HydroModels.step_func

#! define the snow pack reservoir
snow_funcs = [
    SimpleFlux([temp, lday] => [pet],
        exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    SimpleFlux([prcp, temp] => [snowfall, rainfall], [Tmin],
        exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    SimpleFlux([snowpack, temp] => [melt], [Tmax, Df],
        exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
]
snow_dfuncs = [StateFlux([snowfall] => [melt], snowpack)]
snow_ele = HydroBucket(:exphydro_snow, funcs=snow_funcs, dfuncs=snow_dfuncs)

#! define the ET NN and Q NN
et_nn = Lux.Chain(
    Lux.Dense(3 => 16, Lux.tanh),
    Lux.Dense(16 => 16, Lux.leakyrelu),
    Lux.Dense(16 => 1, Lux.leakyrelu),
    name=:etnn
)
q_nn = Lux.Chain(
    Lux.Dense(2 => 16, Lux.tanh),
    Lux.Dense(16 => 16, Lux.leakyrelu),
    Lux.Dense(16 => 1, Lux.leakyrelu),
    name=:qnn
)

#! get init parameters for each NN
et_nn_flux = NeuralFlux([norm_snw, norm_slw, norm_temp] => [log_evap_div_lday], et_nn)
q_nn_flux = NeuralFlux([norm_slw, norm_prcp] => [log_flow], q_nn)

#! define the soil water reservoir
soil_funcs = [
    #* normalize
    StdMeanNormFlux(
        [snowpack, soilwater, prcp, temp] => [norm_snw, norm_slw, norm_prcp, norm_temp],
        [[snowpack_mean, snowpack_std], [soilwater_mean, soilwater_std], [prcp_mean, prcp_std], [temp_mean, temp_std]]
    ),
    et_nn_flux,
    q_nn_flux,
]

state_expr = rainfall + melt - step_func(soilwater) * lday * log_evap_div_lday - step_func(soilwater) * exp(log_flow)
soil_dfuncs = [StateFlux([soilwater, rainfall, melt, lday, log_evap_div_lday, log_flow], soilwater, Num[], expr=state_expr)]
soil_ele = HydroBucket(:m50_soil, funcs=soil_funcs, dfuncs=soil_dfuncs)
soil_ele.ode_func

#! define the Exp-Hydro model
m50_model = HydroModel(:m50, components=[snow_ele, soil_ele]);
#! pretrain each NN
et_nn_input = (norm_snw=snowpack_norm_vec, norm_slw=soilwater_norm_vec, norm_temp=temp_norm_vec)
q_nn_input = (norm_slw=soilwater_norm_vec, norm_prcp=prcp_norm_vec)

et_nn_p = Vector(ComponentVector(LuxCore.initialparameters(StableRNG(42), et_nn)))
q_nn_p = Vector(ComponentVector(LuxCore.initialparameters(StableRNG(42), q_nn)))

input = (prcp=prcp_vec, lday=dayl_vec, temp=temp_vec)
params = (f=0.0167, Smax=1709.46, Qmax=18.47, Df=2.674, Tmax=0.17, Tmin=-2.09)
var_stds = NamedTuple{Tuple([Symbol(nm, :_std) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(stds)
var_means = NamedTuple{Tuple([Symbol(nm, :_mean) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(means)

nn_params = (etnn=et_nn_p, qnn=q_nn_p)
pas = ComponentVector(initstates=(snowpack=0.0, soilwater=1303.00), params=reduce(merge, [params, var_means, var_stds]), nn=nn_params)
timeidx = collect(1:length(prcp_vec))
solver = HydroModels.ODESolver(alg=Tsit5(), reltol=1e-3, abstol=1e-3)
result = m50_model(input, pas, timeidx=timeidx, solver=solver)

# et_nn_p_trained = HydroModels.nn_param_optim(
#     et_nn_flux,
#     input=et_nn_input,
#     target=(log_evap_div_lday=log.(df[ts, :evap] ./ df[ts, :lday]),),
#     init_params=ComponentVector(etnn=et_nn_p),
#     maxiters=100,
# )

# q_nn_p_trained = HydroModels.nn_param_optim(
#     q_nn_flux,
#     input=q_nn_input,
#     target=(log_flow=log.(df[ts, :flow]),),
#     init_params=ComponentVector(qnn=q_nn_p),
#     maxiters=100,
# )

#! set the tunable parameters and constant parameters
#! 当仅优化部分参数如nn_params时就会出错
tunable_pas = ComponentVector(nn=(etnn=collect(ComponentVector(et_nn_p)), qnn=collect(ComponentVector(q_nn_p))))
const_pas = ComponentVector(initstates=(snowpack=0.1, soilwater=1303.00), params=reduce(merge, [params, var_means, var_stds]))
default_model_pas = ComponentArray(merge_recursive(NamedTuple(tunable_pas), NamedTuple(const_pas)))
# new_pas = merge_ca(default_model_pas, tunable_pas)
#! prepare flow
output = (log_flow=qobs_vec,)
#! model calibration
best_pas = HydroModels.param_grad_optim(
    m50_model,
    tunable_pas=tunable_pas,
    const_pas=const_pas,
    input=[input],
    target=[output],
    timeidx=[ts],
    adtype=Optimization.AutoZygote(),
    maxiters=100
)

#! use the optimized parameters for model simulation
# result_opt = m50_model(input, HydroModels.update_ca(default_model_pas, best_pas), timeidx=timeidx, solver=solver)
# pas_list = [0.0, 1303.0042478479704, 0.0167447802633775, 1709.4610152413964, 18.46996175240424, 2.674548847651345, 0.17573919612506747, -2.0929590840638728, 0.8137969540102923]
# pas = ComponentVector(pas_list, getaxes(tunable_pas))

# HydroModels.mse(result_opt.log_flow, qobs_vec)
# 1 - HydroModels.nse(result_opt.log_flow, qobs_vec)

