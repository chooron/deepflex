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
using CairoMakie
using BenchmarkTools
using StableRNGs
using RecursiveArrayTools
include("../../src/LumpedHydro.jl")

# load data
df = DataFrame(CSV.File("temp.csv"));
prcp_vec = df[!, "prcp"]
temp_vec = df[!, "temp"]
dayl_vec = df[!, "lday"]
snowpack_vec = df[!, "snowpack"]
soilwater_vec = df[!, "soilwater"]
qobs_vec = df[!, "flow"]

inputs = [prcp_vec, temp_vec, snowpack_vec, soilwater_vec]
means = mean.(inputs)
stds = std.(inputs)
prcp_norm_vec, temp_norm_vec, snowpack_norm_vec, soilwater_norm_vec = [@.((tmp_vec - mean) / std) for (tmp_vec, mean, std) in zip(inputs, means, stds)]

#! parameters in the Exp-Hydro model
@parameters Tmin = 0.0 [description = "snowfall temperature", unit = "°C"]
@parameters Tmax = 0.0 [description = "snowmelt temperature", unit = "°C"]
@parameters Df = 0.0 [description = "thermal degree-day factor", unit = "mm/(d°C)"]
@parameters Smax = 0.0 [description = "maximum water storage", unit = "mm"]
@parameters f = 0.0 [description = "runoff decline rate", unit = "mm^(-1)"]
@parameters Qmax = 0.0 [description = "maximum subsurface runoff", unit = "mm/d"]
@parameters mean_snowpack = 0.0
@parameters std_snowpack = 0.0
@parameters mean_soilwater = 0.0
@parameters std_soilwater = 0.0
@parameters mean_prcp = 0.0
@parameters std_prcp = 0.0
@parameters mean_temp = 0.0
@parameters std_temp = 0.0
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
snow_dfuncs = [StateFlux([snowfall] => [melt], snowpack, flux_funcs=snow_funcs)]
snow_ele = HydroElement(:exphydro_snow, funcs=snow_funcs, dfuncs=snow_dfuncs)

#! define the ET NN and Q NN
et_nn = Lux.Chain(
    Lux.Dense(3 => 16), # , Lux.tanh
    # Lux.Dense(16 => 16), # , Lux.leakyrelu
    Lux.Dense(16 => 2)#, Lux.leakyrelu
)

q_nn = Lux.Chain(
    Lux.Dense(2 => 16),#, Lux.tanh
    # Lux.Dense(16 => 16),#, Lux.leakyrelu
    Lux.Dense(16 => 1)#, Lux.leakyrelu
)

#! get init parameters for each NN
et_nn_p = LuxCore.initialparameters(StableRNG(42), et_nn)
q_nn_p = LuxCore.initialparameters(StableRNG(42), q_nn)

init_params = Lux.initialparameters(StableRNG(42), et_nn)

chain_params = first(@parameters p[1:length(init_params)] = Vector(ComponentVector(init_params)))
chain_params_type_var = first(@parameters ptype::typeof(typeof(init_params)) = typeof(init_params) [tunable = false])
lazyconvert_params = Symbolics.array_term(convert, chain_params_type_var, chain_params, size=size(chain_params))

@named input = RealInputArray(nin = 3)
exprs = LuxCore.stateless_apply(et_nn, [norm_snw, norm_slw, norm_temp], Symbolics.array_term(convert, ptype, p, size = size(p)))[1]

et_nn_flux = NeuralFlux([norm_snw, norm_slw, norm_temp] => [log_evap_div_lday, log_flow], :etnn => et_nn)
q_nn_flux = NeuralFlux([norm_slw, norm_prcp] => [log_flow], :qnn => q_nn)

et_nn_input = (norm_snw=snowpack_norm_vec, norm_slw=soilwater_norm_vec, norm_temp=temp_norm_vec)
q_nn_input = (norm_slw=soilwater_norm_vec, norm_prcp=prcp_norm_vec)

et_nn_flux_expr = Symbolics.rhss(et_nn_flux.flux_eqs)[1]

temp_func1 = build_function(et_nn_flux_expr, [norm_snw, norm_slw, norm_temp], collect(et_nn_flux.param_info), expression=Val{true})

#! define the soil water reservoir
soil_funcs = [
    #* normalize
    StdMeanNormFlux(
        [snowpack, soilwater, prcp, temp] => [norm_snw, norm_slw, norm_prcp, norm_temp],
        [[mean_snowpack, std_snowpack], [mean_soilwater, std_soilwater], [mean_prcp, std_prcp], [mean_temp, std_temp]]
    ),
    et_nn_flux,
    q_nn_flux,
]

state_expr = rainfall + melt - step_func(soilwater) * lday * log_evap_div_lday - step_func(soilwater) * exp(log_flow)
soil_dfuncs = [StateFlux([soilwater, rainfall, melt, lday, log_evap_div_lday, log_flow], soilwater, Num[], flux_funcs=soil_funcs, state_expr=state_expr)]
soil_ele = HydroElement(:m50_soil, funcs=soil_funcs, dfuncs=soil_dfuncs)

# #! define the Exp-Hydro model
m50_model = HydroUnit(:m50, components=[snow_ele, soil_ele]);
soil_dfuncs[1].inner_func
# #! pretrain each NN
# et_nn_p_trained = LumpedHydro.nn_param_optim(
#     et_nn_flux,
#     input=et_nn_input,
#     target=(log_evap_div_lday=log.(df[!, :evap] ./ df[!, :lday]),),
#     init_params=ComponentVector(etnn=et_nn_p),
#     maxiters=100,
# )

# q_nn_p_trained = LumpedHydro.nn_param_optim(
#     q_nn_flux,
#     input=q_nn_input,
#     target=(log_flow=log.(df[!, :flow]),),
#     init_params=ComponentVector(qnn=q_nn_p),
#     maxiters=100,
# )

# prepare args
# input = (prcp=prcp_vec, lday=dayl_vec, temp=temp_vec)
# params = (f=0.0167, Smax=1709.46, Qmax=18.47, Df=2.674, Tmax=0.17, Tmin=-2.09, etnn=et_nn_p, qnn=q_nn_p)
# mean_params = NamedTuple{Tuple([Symbol(:mean_, nm) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(means)
# std_params = NamedTuple{Tuple([Symbol(:std_, nm) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(stds)
# params = reduce(merge, [params, mean_params, std_params])
# pas = ComponentVector((initstates=(snowpack=0.0, soilwater=1303.00), params=params))
# param_list = collect(pas[:params])
# timeidx = collect(1:length(prcp_vec))
# solver = LumpedHydro.ODESolver(alg=Tsit5(), reltol=1e-3, abstol=1e-3, saveat=timeidx)
# result = m50_model(input, pas, timeidx=timeidx, solver=solver);

# LumpedHydro.get_input_names(soil_ele)
# #! set the tunable parameters and constant parameters
# tunable_pas = ComponentVector(initstates=(snowpack=0.1, soilwater=1303.00),
#     params=(f=0.0167, Smax=1709.46, Qmax=18.47, Df=2.674, Tmax=0.17, Tmin=-2.09))
# const_pas = ComponentVector()
# #! set the tunable parameters boundary
# lower_bounds = [0.01, 100.0, 0.0, 100.0, 10.0, 0.01, 0.0, -3.0]
# upper_bounds = [1500.0, 1500.0, 0.1, 2000.0, 50.0, 5.0, 3.0, 0.0]
# #! prepare flow
# output = (flow=qobs_vec,)
# #! model calibration
# best_pas = LumpedHydro.param_box_optim(
#     exphydro_model,
#     tunable_pas=tunable_pas,
#     const_pas=const_pas,
#     input=input,
#     target=output,
#     target_name=:flow,
#     timeidx=timeidx,
#     lb=lower_bounds,
#     ub=upper_bounds,
#     solve_alg=BBO_adaptive_de_rand_1_bin_radiuslimited(),
#     maxiters=1000,
#     loss_func=LumpedHydro.nse,
#     solver=solver
# )
# #! use the optimized parameters for model simulation
# result_opt = exphydro_model(input, ComponentVector(best_pas; const_pas...), timeidx=timeidx, solver=solver)
# pas_list = [0.0, 1303.0042478479704, 0.0167447802633775, 1709.4610152413964, 18.46996175240424, 2.674548847651345, 0.17573919612506747, -2.0929590840638728, 0.8137969540102923]
# pas = ComponentVector(pas_list, getaxes(tunable_pas))

# 1 - LumpedHydro.nse(result_opt.flow, qobs_vec)
# 1 - LumpedHydro.nse(result.flow, qobs_vec)

