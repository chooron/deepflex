# import lib
using CSV
using Random
using DataFrames
using CairoMakie
using ComponentArrays
using CairoMakie: Axis
using Interpolations

# test gr4j model
include("../../src/DeepFlex.jl")

seed = 42
Random.seed!(42)
P = zeros(100)
E = zeros(100)
P[1:10] = rand(0.0:0.01:10.0, 10)
P[26:30] = rand(0.0:0.01:20.0, 5)
P[41:60] = rand(0.0:0.01:5.0, 20)
P[81:83] = rand(30.0:0.01:50.0, 3)

# build model
s_max, b, a, kf, ks = 10, 5, 0.5, 0.5, 0.5
params = (Smax=s_max, b=b, a=a, kf=kf, ks=ks)
init_states = (soilwater=0.0, fastrouting1=1.0, fastrouting2=1.0,
    fastrouting3=1.0, slowrouting=1.0
)

model = DeepFlex.HyMOD(name=:hymod)

input = (prcp=P, pet=E,time=1:length(E))
output = model(input, params, init_states)

# plot result
fig = Figure(size=(400, 300))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
x = range(1, 100, length=100)
lines!(ax, x, output[:flow], color=:blue)
fig