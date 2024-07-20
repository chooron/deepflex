using CSV
using Lux
using LuxCore
using Random
using DataFrames
using Symbolics
using ComponentArrays
using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using BenchmarkTools
using StableRNGs
using Optimization
using NamedTupleTools
include("../../src/LumpedHydro.jl")

# load data
df = DataFrame(CSV.File("temp.csv"));
ts = collect(1:100)
prcp_vec = df[ts, "prcp"]
temp_vec = df[ts, "temp"]
dayl_vec = df[ts, "lday"]
snowpack_vec = df[ts, "snowpack"]
soilwater_vec = df[ts, "soilwater"]
qobs_vec = df[ts, "flow"]

inputs = [prcp_vec, temp_vec, snowpack_vec, soilwater_vec]
means = mean.(inputs)
stds = std.(inputs)
prcp_norm_vec, temp_norm_vec, snowpack_norm_vec, soilwater_norm_vec =
    [@.((tmp_vec - mean) / std) for (tmp_vec, mean, std) in zip(inputs, means, stds)]

#! parameters in the Exp-Hydro model
@parameters Tmin = 0.0 [description = "snowfall temperature", unit = "°C"]
@parameters Tmax = 0.0 [description = "snowmelt temperature", unit = "°C"]
@parameters Df = 0.0 [description = "thermal degree-day factor", unit = "mm/(d°C)"]
@parameters Smax = 0.0 [description = "maximum water storage", unit = "mm"]
@parameters f = 0.0 [description = "runoff decline rate", unit = "mm^(-1)"]
@parameters Qmax = 0.0 [description = "maximum subsurface runoff", unit = "mm/d"]
@parameters snowpack_norm_param[1:2]
@parameters soilwater_norm_param[1:2]
@parameters prcp_norm_param[1:2]
@parameters temp_norm_param[1:2]

#! hydrological flux in the Exp-Hydro model
@variables prcp(t) = 0.0 [description = "precipitation", unit = "mm"]
@variables temp(t) = 0.0 [description = "precipitation", unit = "°C"]
@variables lday(t) = 0.0 [description = "length of day", unit = "-"]
@variables pet(t) = 0.0 [description = "potential evapotranspiration", unit = "mm"]
@variables snowpack(t) = 0.0 [description = "snow storage", unit = "mm"]
@variables soilwater(t) = 0.0 [description = "catchment water storage", unit = "mm"]
@variables rainfall(t) = 0.0 [description = "rain splitted from precipitation", unit = "mm"]
@variables snowfall(t) = 0.0 [description = "snow splitted from precipitation", unit = "mm"]

@variables melt(t) = 0.0 [description = "melting", unit = "mm"]
@variables log_evap_div_lday(t) = 0.0 [description = "log(evap/lday)"]
@variables log_flow(t) = 0.0 [description = "log discharge"]

@variables norm_snw(t) = 0.0 [description = "discharge", unit = "mm"]
@variables norm_slw(t) = 0.0 [description = "discharge", unit = "mm"]
@variables norm_temp(t) = 0.0 [description = "discharge", unit = "mm"]
@variables norm_prcp(t) = 0.0 [description = "discharge", unit = "mm"]

SimpleFlux = LumpedHydro.SimpleFlux
StdMeanNormFlux = LumpedHydro.StdMeanNormFlux
NeuralFlux = LumpedHydro.NeuralFlux
LagFlux = LumpedHydro.LagFlux
StateFlux = LumpedHydro.StateFlux
HydroElement = LumpedHydro.HydroElement
HydroUnit = LumpedHydro.HydroUnit
step_func = LumpedHydro.step_func

#! define the snow pack reservoir
snow_funcs = [
    SimpleFlux([temp, lday] => [pet],
        flux_exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    SimpleFlux([prcp, temp] => [snowfall, rainfall], [Tmin],
        flux_exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    SimpleFlux([snowpack, temp] => [melt], [Tmax, Df],
        flux_exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
]
snow_dfuncs = [StateFlux([snowfall] => [melt], snowpack)]
snow_ele = HydroElement(:exphydro_snow, funcs=snow_funcs, dfuncs=snow_dfuncs)

#! define the ET NN and Q NN
et_nn = Lux.Chain(
    Lux.Dense(3 => 16, Lux.tanh),
    Lux.Dense(16 => 16, Lux.leakyrelu),
    Lux.Dense(16 => 2, Lux.leakyrelu)
)
q_nn = Lux.Chain(
    Lux.Dense(2 => 16, Lux.tanh),
    Lux.Dense(16 => 16, Lux.leakyrelu),
    Lux.Dense(16 => 1, Lux.leakyrelu)
)

#! get init parameters for each NN
et_nn_flux = NeuralFlux([norm_snw, norm_slw, norm_temp] => [log_evap_div_lday, log_flow], :etnn => et_nn)
q_nn_flux = NeuralFlux([norm_slw, norm_prcp] => [log_flow], :qnn => q_nn)

#! define the soil water reservoir
soil_funcs = [
    #* normalize
    StdMeanNormFlux(
        [snowpack, soilwater, prcp, temp] => [norm_snw, norm_slw, norm_prcp, norm_temp],
        [snowpack_norm_param, soilwater_norm_param, prcp_norm_param, temp_norm_param]
    ),
    NeuralFlux([norm_snw, norm_slw, norm_temp] => [log_evap_div_lday, log_flow], :etnn => et_nn),
    NeuralFlux([norm_slw, norm_prcp] => [log_flow], :qnn => q_nn),
]

state_expr = rainfall + melt - step_func(soilwater) * lday * log_evap_div_lday - step_func(soilwater) * exp(log_flow)
soil_dfuncs = [StateFlux([soilwater, rainfall, melt, lday, log_evap_div_lday, log_flow], soilwater, Num[], state_expr=state_expr)]
soil_ele = HydroElement(:m50_soil, funcs=soil_funcs, dfuncs=soil_dfuncs)

#! define the Exp-Hydro model
m50_model = HydroUnit(:m50, components=[snow_ele, soil_ele]);


#! pretrain each NN
et_nn_input = (norm_snw=snowpack_norm_vec, norm_slw=soilwater_norm_vec, norm_temp=temp_norm_vec)
q_nn_input = (norm_slw=soilwater_norm_vec, norm_prcp=prcp_norm_vec)

# et_nn_p_trained = LumpedHydro.nn_param_optim(
#     et_nn_flux,
#     input=et_nn_input,
#     target=(log_evap_div_lday=log.(df[ts, :evap] ./ df[ts, :lday]),),
#     init_params=ComponentVector(etnn=et_nn_p),
#     maxiters=100,
# )

# q_nn_p_trained = LumpedHydro.nn_param_optim(
#     q_nn_flux,
#     input=q_nn_input,
#     target=(log_flow=log.(df[ts, :flow]),),
#     init_params=ComponentVector(qnn=q_nn_p),
#     maxiters=100,
# )

# prepare args
et_nn_p = LuxCore.initialparameters(StableRNG(42), et_nn)
q_nn_p = LuxCore.initialparameters(StableRNG(42), q_nn)

input = (prcp=prcp_vec, lday=dayl_vec, temp=temp_vec)
params = (f=0.0167, Smax=1709.46, Qmax=18.47, Df=2.674, Tmax=0.17, Tmin=-2.09)
norm_params = NamedTuple{Tuple([Symbol(nm, :_norm_param) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(
    [[mean, std] for (mean, std) in zip(means, stds)]
)
nn_params = (etnn=et_nn_p, qnn=q_nn_p)
pas = ComponentVector((initstates=(snowpack=0.0, soilwater=1303.00), params=reduce(merge, [params, norm_params, nn_params])))
timeidx = collect(1:length(prcp_vec))
solver = LumpedHydro.ODESolver(alg=Tsit5(), reltol=1e-3, abstol=1e-3, saveat=timeidx)
# result = m50_model(input, pas, timeidx=timeidx, solver=solver)

#! set the tunable parameters and constant parameters
#! 当仅优化部分参数如nn_params时就会出错
tunable_pas = ComponentVector(params=(etnn=collect(ComponentVector(et_nn_p)), qnn=collect(ComponentVector(q_nn_p))))
const_pas = ComponentVector(initstates=(snowpack=0.1, soilwater=1303.00), params=reduce(merge, [params, norm_params]))
default_model_pas = ComponentArray(merge_recursive(NamedTuple(tunable_pas), NamedTuple(const_pas)))
# new_pas = merge_ca(default_model_pas, tunable_pas)
#! prepare flow
output = (log_flow=qobs_vec,)
#! model calibration
best_pas = LumpedHydro.param_grad_optim(
    m50_model,
    tunable_pas=tunable_pas,
    const_pas=const_pas,
    input=input,
    target=output,
    timeidx=ts,
    adtype=Optimization.AutoZygote()
)

#! use the optimized parameters for model simulation
result_opt = m50_model(input, LumpedHydro.merge_ca(default_model_pas, best_pas), timeidx=timeidx, solver=solver)
pas_list = [0.0, 1303.0042478479704, 0.0167447802633775, 1709.4610152413964, 18.46996175240424, 2.674548847651345, 0.17573919612506747, -2.0929590840638728, 0.8137969540102923]
pas = ComponentVector(pas_list, getaxes(tunable_pas))

LumpedHydro.mse(result_opt.log_flow, qobs_vec)
1 - LumpedHydro.nse(result.flow, qobs_vec)

