# import lib
using CSV
using DataFrames
using ComponentArrays
using Zygote, Lux
using BenchmarkTools
using Random
using ModelingToolkit

@variables A[1:10, 1:10]

# test exphydro model
# include("../../src/LumpedHydro.jl")

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

lux_model = Lux.Chain(Lux.Dense(3, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 1, leakyrelu),name=:model)
rng = MersenneTwister()
Random.seed!(rng, 42)
ps,st = Lux.setup(rng, lux_model)
ps = Lux.initialparameters(rng,lux_model)
ps
# todo 后续需要将lux的参数嵌入至mtk中
# initial_params(lux_model)
lux_func = LumpedHydro.LuxNNFlux([:prcp, :temp, :lday], [:flow], lux_model=lux_model)
# @btime LumpedHydro.nn_param_optim(lux_func, input=inputs, target=outputs, init_params=lux_func.init_params)
using DiffEqFlux

ann = FastChain(FastDense(1,32,tanh), FastDense(32,32,tanh), FastDense(32,1))