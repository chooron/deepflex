using CSV
using Random
using DataFrames
using ComponentArrays
using OrdinaryDiffEq
using ModelingToolkit
using OptimizationBBO
using Optimization
using BenchmarkTools
using DataInterpolations
using Plots
using HydroErrors
include("../../src/HydroModels.jl")

#! parameters in the Exp-Hydro model
@parameters Tmin = 0.0 [description = "snowfall temperature", unit = "°C"]
@parameters Tmax = 0.0 [description = "snowmelt temperature", unit = "°C"]
@parameters Df = 0.0 [description = "thermal degree-day factor", unit = "mm/(d°C)"]
@parameters Smax = 0.0 [description = "maximum water storage", unit = "mm"]
@parameters f = 0.0 [description = "runoff decline rate", unit = "mm^(-1)"]
@parameters Qmax = 0.0 [description = "maximum subsurface runoff", unit = "mm/d"]
#! hydrological flux in the Exp-Hydro model
@variables prcp = 0.0 [description = "precipitation", unit = "mm"]
@variables temp = 0.0 [description = "precipitation", unit = "°C"]
@variables lday = 0.0 [description = "length of day", unit = "-"]
@variables pet = 0.0 [description = "potential evapotranspiration", unit = "mm"]
@variables snowpack = 0.0 [description = "snow storage", unit = "mm"]
@variables soilwater = 0.0 [description = "catchment water storage", unit = "mm"]
@variables rainfall = 0.0 [description = "rain splitted from precipitation", unit = "mm"]
@variables snowfall = 0.0 [description = "snow splitted from precipitation", unit = "mm"]
@variables evap = 0.0 [description = "evapotranspiration", unit = "mm"]
@variables melt = 0.0 [description = "melting", unit = "mm"]
@variables baseflow = 0.0 [description = "base discharge", unit = "mm"]
@variables surfaceflow = 0.0 [description = "surface discharge", unit = "mm"]
@variables flow = 0.0 [description = "discharge", unit = "mm"]
HydroFlux = HydroModels.HydroFlux
StateFlux = HydroModels.StateFlux
HydroBucket = HydroModels.HydroBucket
HydroModel = HydroModels.HydroModel
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5
# HydroFlux([temp, lday] => [pet],
# exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
#! define the snow pack reservoir

snow_funcs = [
    HydroFlux([temp, lday] => [pet],
        exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin],
        exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    HydroFlux([snowpack, temp] => [melt], [Tmax, Df],
        exprs=[step_func(temp - Tmax) * min(snowpack, Df * (temp - Tmax))]),
]
snow_dfuncs = [StateFlux([snowfall] => [melt], snowpack)]

snow_ele = HydroBucket(name=:exphydro_snow, funcs=snow_funcs, dfuncs=snow_dfuncs)
snow_ele.flux_func
#! define the soil water reservoir
soil_funcs = [
    HydroFlux([soilwater, pet] => [evap], [Smax],
        exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
    HydroFlux([soilwater] => [baseflow], [Smax, Qmax, f],
        exprs=[step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))]),
    HydroFlux([soilwater] => [surfaceflow], [Smax],
        exprs=[max(0.0, soilwater - Smax)]),
    HydroFlux([baseflow, surfaceflow] => [flow],
        exprs=[baseflow + surfaceflow]),
]

soil_funcs = [
    HydroFlux([soilwater, pet] => [evap], [Smax],
        exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
    HydroFlux([soilwater] => [flow], [Smax, Qmax, f],
        exprs=[Qmax * exp(-f * (Smax - soilwater)) +
               max(0.0, soilwater - Smax)]),
]
soil_dfuncs = [StateFlux([rainfall, melt] => [evap, flow], soilwater)]

soil_ele = HydroBucket(
    name=:exphydro_soil,
    funcs=soil_funcs,
    dfuncs=soil_dfuncs
)

#! define the Exp-Hydro model
exphydro_model = HydroModel(
    name=:exphydro,
    components=[snow_ele, soil_ele]
)

# load data
df = DataFrame(CSV.File("data/exphydro/01013500.csv"));
prcp_vec = df[!, "prcp(mm/day)"]
temp_vec = df[!, "tmean(C)"]
dayl_vec = df[!, "dayl(day)"]
qobs_vec = df[!, "flow(mm)"]

# prepare args
input = (prcp=prcp_vec, lday=dayl_vec, temp=temp_vec)
pas = ComponentVector((initstates=(snowpack=0.0, soilwater=1303.00),
    params=(f=0.0167, Smax=1709.46, Qmax=18.47, Df=2.674, Tmax=0.17, Tmin=-2.09)))
timeidx = collect(1:length(prcp_vec))
solver = HydroModels.ODESolver(alg=Tsit5(), reltol=1e-3, abstol=1e-3, saveat=timeidx)
config = (solver=solver, timeidx=timeidx)
@btime result = exphydro_model(input, pas, config=config, convert_to_ntp=false); #  29.317 ms (356101 allocations: 28.18 MiB)

# plot(result.flow)
# plot!(qobs_vec)
# #! set the tunable parameters and constant parameters
# tunable_pas = ComponentVector(params=(f=0.0167, Smax=1709.46, Qmax=18.47, Df=2.674, Tmax=0.17, Tmin=-2.09))
# const_pas = ComponentVector(initstates=(snowpack=0.1, soilwater=1303.00))
# #! set the tunable parameters boundary
# lower_bounds = [0.0, 100.0, 10.0, 0.01, 0.0, -3.0]
# upper_bounds = [0.1, 2000.0, 50.0, 5.0, 3.0, 0.0]
# #! prepare flow
# output = (flow=qobs_vec,)
# #! model calibration
# # best_pas = HydroModels.param_box_optim(
# #     exphydro_model,
# #     tunable_pas=tunable_pas,
# #     const_pas=const_pas,
# #     input=[input],
# #     target=[output],
# #     lb=lower_bounds,
# #     ub=upper_bounds,
# #     solve_alg=BBO_adaptive_de_rand_1_bin_radiuslimited(),
# #     maxiters=1000,
# #     loss_func=HydroErrors.mse,
# # )
# #! use the optimized parameters for model simulation
# best_pas = HydroModels.param_grad_optim(
#     exphydro_model,
#     tunable_pas=tunable_pas,
#     const_pas=const_pas,
#     input=[input],
#     target=[output],
#     config=[(solver=HydroModels.DiscreteSolver(),)],
#     maxiters=100,
#     loss_func=HydroErrors.mse,
#     adtype=Optimization.AutoForwardDiff(),
# )

# total_pas = ComponentVector(params=(f=0.0167, Smax=1709.46, Qmax=18.47, Df=2.674, Tmax=0.17, Tmin=-2.09), initstates=(snowpack=0.1, soilwater=1303.00))
# update_pas = HydroModels.update_ca(total_pas, ComponentVector(best_pas, getaxes(tunable_pas)))
# result_opt = exphydro_model(input, update_pas, timeidx=timeidx, solver=solver)
# result_opt_df = DataFrame(result_opt)
# result_opt_df.qobs = qobs_vec
# # CSV.write("temp.csv", result_opt_df)
# 1 - HydroErr.nse(result_opt.flow, qobs_vec)
# 1 - HydroErr.nse(result.flow, qobs_vec)

