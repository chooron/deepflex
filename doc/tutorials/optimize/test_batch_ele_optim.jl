# import lib
using CSV
using DataFrames
using ComponentArrays
using OptimizationOptimisers
using BenchmarkTools
using NamedTupleTools
using Optimization
using HydroErrors

include("../../../src/HydroModels.jl")

# test exphydro model

# predefine the parameters
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

tunable_pas = ComponentVector(params=(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin))
const_pas = ComponentVector(initstates=(snowpack=0.0, soilwater=1300.0))

init_pas = ComponentVector(params=(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin), initstates=(snowpack=0.0, soilwater=1300.0))
model = HydroModels.ExpHydro.SurfaceStorage(name=:sf)

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
input_matrix = Matrix(reduce(hcat, [input[k] for k in keys(input)])')
run_kwargs = (convert_to_ntp=true,)

best_pas = HydroModels.param_batch_optim(
    model,
    tunable_pas=tunable_pas,
    const_pas=const_pas,
    input=repeat([input], 10),
    target=repeat([output], 10),
    config=fill((solver=HydroModels.ODESolver(),timeidx=timeidx), 10),
    adtype=Optimization.AutoZygote(),
    maxiters=10,
    loss_func=HydroErrors.mse,
    run_kwargs=(convert_to_ntp=true,)
)