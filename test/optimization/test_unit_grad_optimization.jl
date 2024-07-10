# import lib
using CSV
using DataFrames
# using CairoMakie
using ComponentArrays
using StructArrays
using OptimizationOptimisers
using BenchmarkTools
using NamedTupleTools
using Optimization
# using LumpedHydro

# test exphydro model
include("../../src/LumpedHydro.jl")

# predefine the parameters
# init_parameter = [0.0, 100.0, 0.01, 20, 1.0, 1.0, -1.0]
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

tunable_pas = ComponentVector(params=ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin))
const_pas = ComponentVector(initstates=ComponentVector(snowwater=0.0, soilwater=1300.0))

model = LumpedHydro.ExpHydro.Unit(name=:exphydro)

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
output = (flow=flow_vec,)

best_pas = LumpedHydro.param_grad_optim(
    model,
    tunable_pas=tunable_pas,
    const_pas=const_pas,
    input=input,
    target=output,
    timeidx=ts,
    adtype=Optimization.AutoZygote()
)