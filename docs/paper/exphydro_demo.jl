using CSV
using Random
using DataFrames
using ComponentArrays
using OrdinaryDiffEq
using ModelingToolkit
using CairoMakie
using OptimizationBBO
using BenchmarkTools
using DataInterpolations
include("../../src/LumpedHydro.jl")

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
@variables baseflow = 0.0 [description = "discharge", unit = "mm"]
@variables surfaceflow = 0.0 [description = "discharge", unit = "mm"]
@variables flow = 0.0 [description = "discharge", unit = "mm"]
SimpleFlux = LumpedHydro.SimpleFlux
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

#! define the soil water reservoir
soil_funcs = [
    SimpleFlux([soilwater, pet] => [evap], [Smax],
        flux_exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
    SimpleFlux([soilwater] => [baseflow], [Smax, Qmax, f],
        flux_exprs=[step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))]),
    SimpleFlux([soilwater] => [surfaceflow], [Smax],
        flux_exprs=[max(0.0, soilwater - Smax)]),
    SimpleFlux([baseflow, surfaceflow] => [flow],
        flux_exprs=[baseflow + surfaceflow]),
]
soil_dfuncs = [StateFlux([rainfall, melt] => [evap, flow], soilwater)]
soil_ele = HydroElement(:exphydro_soil, funcs=soil_funcs, dfuncs=soil_dfuncs)

#! define the Exp-Hydro model
exphydro_model = HydroUnit(:exphydro, components=[snow_ele, soil_ele])

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
solver = LumpedHydro.ODESolver(alg=Tsit5(), reltol=1e-3, abstol=1e-3, saveat=timeidx)
result = exphydro_model(input, pas, timeidx=timeidx, solver=solver);

@info 1 - LumpedHydro.nse(result.flow, qobs_vec)

fig = Figure(size=(400, 300))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
lines!(ax, timeidx, result.flow, color=:red)
lines!(ax, timeidx, qobs_vec, color=:blue)
fig

#! set the tunable parameters and constant parameters
tunable_pas = ComponentVector(params=(f=0.0167, Smax=1709.46, Qmax=18.47, Df=2.674, Tmax=0.17, Tmin=-2.09))
const_pas = ComponentVector(initstates=(snowpack=0.1, soilwater=1303.00))
#! set the tunable parameters boundary
lower_bounds = [0.0, 100.0, 10.0, 0.01, 0.0, -3.0]
upper_bounds = [0.1, 2000.0, 50.0, 5.0, 3.0, 0.0]
#! prepare flow
output = (flow=qobs_vec,)
#! model calibration
best_pas = LumpedHydro.param_box_optim(
    exphydro_model,
    tunable_pas=tunable_pas,
    const_pas=const_pas,
    input=input,
    target=output,
    timeidx=timeidx,
    lb=lower_bounds,
    ub=upper_bounds,
    solve_alg=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    maxiters=1000,
    loss_func=LumpedHydro.nse,
    solver=solver
)
# #! use the optimized parameters for model simulation
# result_opt = exphydro_model(input, ComponentVector(best_pas; const_pas...), timeidx=timeidx, solver=solver)
# pas_list = [0.0, 1303.0042478479704, 0.0167447802633775, 1709.4610152413964, 18.46996175240424, 2.674548847651345, 0.17573919612506747, -2.0929590840638728, 0.8137969540102923]
# pas = ComponentVector(pas_list, getaxes(tunable_pas))
# result_opt_df = DataFrame(result_opt)
# result_opt_df.qobs = qobs_vec
# CSV.write("temp.csv", result_opt_df)
# 1 - LumpedHydro.nse(result_opt.flow, qobs_vec)
# 1 - LumpedHydro.nse(result.flow, qobs_vec)

