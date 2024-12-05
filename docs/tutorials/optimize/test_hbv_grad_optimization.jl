# import lib
using CSV
using DataFrames
# using CairoMakie
using ComponentArrays
using OptimizationOptimisers
using BenchmarkTools
using NamedTupleTools
using Optimization
using HydroModels
# test exphydro model
# include("../../src/HydroModels.jl")
# 5.0.0
# predefine the parameters
tt, dd, FC, Beta, PWP, L = 0.0, 4.25, 177.1, 2.35, 105.89, 4.87
k0, k1, k2, kp = 0.05, 0.03, 0.02, 0.05
params = ComponentVector(tt=tt, dd=dd, FC=FC, Beta=Beta, PWP=PWP, L=L, k0=k0, k1=k1, k2=k2, kp=kp)
initstates = ComponentVector(snowwater=0.0, soilwater=100.0, s1=3, s2=10)
pas = ComponentVector((vertical=(params=params, initstates=initstates)))

tunable_pas = ComponentVector(params=params)
const_pas = ComponentVector(initstates=initstates)

unit = HydroModels.HBV_EDU.Model(name=:hbv)

# load data
file_path = "data/hbv_edu/hbv_sample.csv"
data = CSV.File(file_path);
df = DataFrame(data);

input = (prcp=df[!,:prec], pet=df[!,:pet], temp=df[!,:temp])
target = (flow=df[!,:qsim],)
timeidx = collect(1:length(df[!,:qsim]))

best_pas = HydroModels.param_grad_optim(
    unit,
    tunable_pas=tunable_pas,
    const_pas=const_pas,
    input=[input],
    target=[target],
    timeidx=[timeidx],
    adtype=Optimization.AutoZygote()
)
