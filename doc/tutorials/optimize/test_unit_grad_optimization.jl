# import lib
using CSV
using DataFrames
using ComponentArrays
using OptimizationOptimisers
using BenchmarkTools
using NamedTupleTools
using Optimization

# test exphydro model
include("../../../src/HydroModels.jl")

# predefine the parameters
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

tunable_pas = ComponentVector(params=ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin))
const_pas = ComponentVector(initstates=ComponentVector(snowpack=0.0, soilwater=1300.0))

model = HydroModels.ExpHydro.Model(name=:exphydro)

# load data
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:100)
lday_vec = df[ts, "dayl(day)"]
prcp_vec = df[ts, "prcp(mm/day)"]
temp_vec = df[ts, "tmean(C)"]
flow_vec = df[ts, "flow(mm)"]

# parameters optimization
input = (prcp=prcp_vec, lday=lday_vec, temp=temp_vec)

input_matrix = Matrix(reduce(hcat, collect(input))')
output = (flow=flow_vec,)
pas = ComponentVector(tunable_pas; const_pas...)
result = model(input_matrix, pas, ts, kwargs=run_kwargs)


best_pas = HydroModels.param_grad_optim(
    model,
    tunable_pas=tunable_pas,
    const_pas=const_pas,
    input=repeat([input], 10),
    target=repeat([output], 10),
    timeidx=repeat([ts], 10),
    adtype=Optimization.AutoZygote(),
    maxiters=10,
    loss_func=HydroErrors.mse,
)