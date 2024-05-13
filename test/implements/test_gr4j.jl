# import lib
# tested in https://github.com/ckrapu/gr4j_theano
using CSV
using Random
using DataFrames
using CairoMakie
using CairoMakie: Axis
using BenchmarkTools
using ModelingToolkit
using OrdinaryDiffEq
using ComponentArrays

# test gr4j model
include("../../src/LumpedHydro.jl")

# load data
file_path = "data/gr4j/sample.csv"
data = CSV.File(file_path);
df = DataFrame(data);
# time = 1:10000
prcp_vec = df[!, "P"]
et_vec = df[!, "ET"]
# flow_vec = df[time, "flow(mm)"]

# build model

unit_params = (x1=320.11, x2=2.42, x3=69.63, x4=1.39, ω=3.5, γ=5.0)
route_params = (x4=x4,)
init_states = (soilwater=0.6 * 320.11, routingstore=0.70 * 69.63)
pas = ComponentVector(gr4j=(params=unit_params, initstates=init_states, weight=1.0))

model = LumpedHydro.GR4J.Node(name=:gr4j, step=true, mtk=false)
solver = LumpedHydro.ODESolver(alg=Tsit5())
input = (gr4j=(prcp=prcp_vec, pet=et_vec, time=1:1:length(prcp_vec)),)
result = model(input, pas, solver=solver);
result_df = DataFrame(result)

# plot result
fig = Figure(size=(400, 300))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
x = range(1, 100, length=100)
lines!(ax, 1:1:length(prcp_vec), result_df[!, :flow], color=:blue)
df[!, "my_pred"] = result_df[!, :flow]
lines!(ax, 1:1:length(prcp_vec), df[!, "obs_Q"], color=:red)
lines!(ax, 1:1:length(prcp_vec), df[!, "modeled_Q"], color=:green)
fig