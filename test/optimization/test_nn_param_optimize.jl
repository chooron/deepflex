# import lib
using CSV
using DataFrames
using CairoMakie
using ComponentArrays
using Zygote, Lux
using BenchmarkTools

# test exphydro model
include("../../src/DeepFlex.jl")

# load data
file_path = "data/camels/01013500.csv"

data = CSV.File(file_path);
df = DataFrame(data);
lday_vec = df[1:1000, "dayl(day)"]
prcp_vec = df[1:1000, "prcp(mm/day)"]
temp_vec = df[1:1000, "tmean(C)"]
flow_vec = df[1:1000, "flow(mm)"]

inputs = (prcp=prcp_vec, temp=temp_vec, lday=lday_vec)
outputs = (flow=flow_vec,)

lux_model = Lux.Chain(Lux.Dense(3, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 1, leakyrelu))
lux_func = DeepFlex.LuxNNFlux([:prcp, :temp, :lday], [:flow], lux_model=lux_model)
@btime DeepFlex.nn_param_optim(lux_func, input=inputs, target=outputs, init_params=lux_func.init_params)