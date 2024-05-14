# import lib
using CSV
using DataFrames
# using CairoMakie
using ComponentArrays
using OptimizationOptimisers
using BenchmarkTools
using NamedTupleTools
using Optimization

# # test exphydro model
include("../../src/LumpedHydro.jl")

# predefine the parameters
# init_parameter = [0.0, 100.0, 0.01, 20, 1.0, 1.0, -1.0]
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

tunable_pas = ComponentVector(params=(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin))
const_pas = ComponentVector(initstates=(snowwater=0.0, soilwater=1300.0))

model = LumpedHydro.ExpHydro.Surface(name=:sf, mtk=true)

# load data
file_path = "data/m50/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
lday_vec = df[1:1000, "Lday"]
prcp_vec = df[1:1000, "Prcp"]
temp_vec = df[1:1000, "Temp"]
flow_vec = df[1:1000, "SnowWater"]

# parameters optimization
input = (prcp=prcp_vec, lday=lday_vec, temp=temp_vec, snowwater=flow_vec, time=1:1:length(lday_vec))
output = (snowwater=flow_vec,)

best_pas = LumpedHydro.param_grad_optim(
    model,
    tunable_pas=tunable_pas,
    const_pas=const_pas,
    input=input,
    target=output,
    adtype=Optimization.AutoForwardDiff(),
    target_name=:snowwater
)


# total_params = DeepFlex.merge_ca(best_pas, const_pas)[:param]
# result = model(input, total_params)