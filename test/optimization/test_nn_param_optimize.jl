# import lib
using CSV
using DataFrames
using CairoMakie
using ComponentArrays
using Zygote, Lux, LuxCUDA

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

inputs = ComponentVector(prcp=prcp_vec, temp=temp_vec, lday=lday_vec)
outputs = ComponentVector(flow=flow_vec)

q_nn = Lux.Chain(Lux.Dense(3, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 1, leakyrelu))
lux_func = DeepFlex.LuxNNFlux([:prcp, :temp, :lday], [:flow], lux_model=q_nn)

opt_params = DeepFlex.nn_params_optimize!(lux_func, input=inputs, output=outputs)
# output = DeepFlex.get_output(ele, input=inputs)