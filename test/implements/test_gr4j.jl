# import lib
using CSV
using Random
using DataFrames
using CairoMakie
using CairoMakie: Axis
using BenchmarkTools
using ModelingToolkit
using OrdinaryDiffEq

# test gr4j model
include("../../src/DeepFlex.jl")

seed = 42
Random.seed!(42)
P = ones(100)
E = zeros(100)
P[1:10] = rand(0.0:0.01:10.0, 10)
P[26:30] = rand(0.0:0.01:20.0, 5)
P[41:60] = rand(0.0:0.01:5.0, 20)
P[81:83] = rand(30.0:0.01:50.0, 3)

# build model
x1, x2, x3, x4 = 50.0, 0.1, 20.0, 3.5
unit_params = (x1=x1, x2=x2, x3=x3, x4=x4, ω=3.5, γ=5.0)
route_params = (x4=x4,)
init_states = (soilwater=10.0, routingstore=10.0)
pas = ComponentVector(gr4j=(params=unit_params, initstates=init_states, weight=1.0))

model = DeepFlex.GR4J.Node(name=:gr4j)
solver = DeepFlex.ODESolver(alg=Tsit5())
input = (gr4j=(prcp=P, pet=E, time=1:1:length(P)),)
result = model(input, pas, solver=solver);
result_df = DataFrame(result)

# plot result
fig = Figure(size=(400, 300))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
x = range(1, 100, length=100)
lines!(ax, x, result_df[!, :flow], color=:blue)
fig