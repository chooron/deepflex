# import lib
using CSV
using DataFrames
using CairoMakie
using Zygote, Lux, LuxCUDA
using OrdinaryDiffEq
using DiffEqFlux

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

input_layer = Lux.Dense(3, 10, relu)
nn_layer = Lux.Chain(
    Lux.Dense(10, 20, relu),
    Lux.Dense(20, 20, relu),
    Lux.Dense(20, 10, relu)
)
node_layer = NeuralODE(nn_layer,
    (0.0f0, 1.0f0),
    Tsit5();
    save_everystep=false,
    saveat=
    reltol=1e-3, abstol=1e-3,
    save_start=false)
output_layer = Lux.Dense(10, 1)


model_no_ode = Lux.Chain(; input_layer, nn_layer, output_layer)
# ele = DeepFlex.LuxElement(model_no_ode, device=dev_gpu)

# res = DeepFlex.node_params_optimize(ele, input=inputs, output=output)
# # output = DeepFlex.get_output(ele, input=inputs)


