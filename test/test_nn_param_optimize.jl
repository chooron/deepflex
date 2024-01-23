# import lib
using CSV
using DataFrames
using CairoMakie
using Zygote, Lux, LuxCUDA

# test exphydro model
include("../src/DeepFlex.jl")

# load data
file_path = "data/camels/01013500.csv"

data = CSV.File(file_path);
df = DataFrame(data);
lday_vec = df[1:1000, "dayl(day)"]
prcp_vec = df[1:1000, "prcp(mm/day)"]
temp_vec = df[1:1000, "tmean(C)"]
flow_vec = df[1:1000, "flow(mm)"]

dev_gpu = Lux.gpu_device()
inputs = Dict(:Prcp => prcp_vec, :Lday => lday_vec, :Temp => temp_vec)
output = Dict(:Q => flow_vec)
ele = DeepFlex.LinearNNElement(3, 1, 64, dev_gpu)
DeepFlex.nn_params_optimize!(ele, input=inputs, output=output)
output = DeepFlex.get_output(ele, input=inputs)