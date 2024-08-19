# import lib
using CSV
using DataFrames
using ComponentArrays
using OptimizationOptimisers
using BenchmarkTools
using NamedTupleTools
using Optimization
using HydroEquations
using ModelingToolkit

include("../../../src/HydroModels.jl")

@variables temp lday prcp snowpack
@variables pet snowfall melt rainfall infiltration
@parameters Tmin Tmax Df

# test exphydro model
#* test state flux
snow_funcs = [
    HydroModels.SimpleFlux([temp, lday] => [pet],
        flux_exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    HydroModels.SimpleFlux([prcp, temp] => [snowfall, rainfall], [Tmin],
        flux_exprs=[HydroModels.step_func(Tmin - temp) * prcp, HydroModels.step_func(temp - Tmin) * prcp]),
    HydroModels.SimpleFlux([snowpack, temp] => [melt], [Tmax, Df],
        flux_exprs=[HydroModels.step_func(temp - Tmax) * HydroModels.step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
]
snow_dfuncs = [HydroModels.StateFlux([snowfall] => [melt], snowpack)]

model = HydroModels.HydroElement(:sf, funcs=snow_funcs, dfuncs=snow_dfuncs)

# predefine the parameters
tunable_pas = ComponentVector(params=(f=0.01674478, Smax=1709.461015, Qmax=18.46996175, Df=2.674548848, Tmax=0.175739196, Tmin=-2.092959084))
const_pas = ComponentVector(initstates=(snowpack=0.0, soilwater=1300.0))
init_pas = ComponentVector(
    params=(f=0.01674478, Smax=1709.461015, Qmax=18.46996175, Df=2.674548848, Tmax=0.175739196, Tmin=-2.092959084),
    initstates=(snowpack=0.0, soilwater=1300.0)
)

# load data
file_path = "data/m50/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
timeidx = collect(1:100)
lday_vec = df[timeidx, "Lday"]
prcp_vec = df[timeidx, "Prcp"]
temp_vec = df[timeidx, "Temp"]
flow_vec = df[timeidx, "SnowWater"]

# parameters optimization
input = (prcp=prcp_vec, lday=lday_vec, temp=temp_vec,)
output = (melt=flow_vec,)

best_pas = HydroModels.param_grad_optim(
    model,
    tunable_pas=tunable_pas,
    const_pas=const_pas,
    input=repeat([input], 10),
    target=repeat([output], 10),
    timeidx=repeat([timeidx], 10),
    adtype=Optimization.AutoZygote(),
    maxiters=100,
    loss_func=HydroEquations.mse
)
# ComponentArray(merge_recursive(NamedTuple(tunable_pas), NamedTuple(const_pas)))