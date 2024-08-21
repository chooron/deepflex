# import lib
using CSV
using DataFrames
using ComponentArrays
using OptimizationBBO
using BenchmarkTools
using NamedTupleTools
using HydroErrors
# using HydroModels

include("../../../src/HydroModels.jl")
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

tunable_pas = ComponentVector(params=ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin),)
const_pas = ComponentVector((initstates=ComponentVector(snowpack=0.0, soilwater=1300.0), weight=1.0))

params_axes = getaxes(tunable_pas)

lb_list = [0.0, 100.0, 10.0, 0.0, 0.0, -3.0]
ub_list = [0.1, 2000.0, 50.0, 5.0, 3.0, 0.0]

ts = collect(1:10000)

tunable_param_lb = ComponentVector(lb_list, getaxes(tunable_pas))
tunable_param_ub = ComponentVector(ub_list, getaxes(tunable_pas))

model = HydroModels.ExpHydro.Model(name=:exphydro)

# load data
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
lday_vec = df[ts, "dayl(day)"]
prcp_vec = df[ts, "prcp(mm/day)"]
temp_vec = df[ts, "tmean(C)"]
flow_vec = df[ts, "flow(mm)"]

# parameters optimization
input = (prcp=prcp_vec, lday=lday_vec, temp=temp_vec)
output = (flow=flow_vec,)

best_pas = HydroModels.param_box_optim(
    model,
    tunable_pas=tunable_pas,
    const_pas=const_pas,
    input=[input],
    target=[output],
    timeidx=[ts],
    lb=tunable_param_lb,
    ub=tunable_param_ub,
    solve_alg=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    maxiters=100,
    loss_func=HydroErr.mse,
)