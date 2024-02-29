# import lib
using CSV
using Random
using DataFrames
using CairoMakie
using ComponentArrays
using CairoMakie: Axis
using Interpolations

# test gr4j model
include("../src/DeepFlex.jl")

seed = 42
Random.seed!(42)
P = zeros(100)
E = zeros(100)
T = zeros(100)
P[1:10] = rand(0.0:0.01:10.0, 10)
P[26:30] = rand(0.0:0.01:20.0, 5)
P[41:60] = rand(0.0:0.01:5.0, 20)
P[81:83] = rand(30.0:0.01:50.0, 3)

T[20:30] .= 1
T[50:70] .= 10
T[70:90] .= 5

# build model
tt, tti, cfr, cfmax, ttm, whc = 0, 4, 0.5, 10, 0, 0.5
cflux, fc, lp, k0, k1, α, β, c = 2.0, 500, 0.5, 0.5, 0.5, 2, 5, 10
maxbas = 5

parameters = ComponentVector(
    tt=tt, tti=tti, cfr=cfr, cfmax=cfmax, ttm=ttm, whc=whc,
    cflux=cflux, fc=fc, lp=lp, k0=k0, k1=k1, α=α, β=β, c=c, maxbas=maxbas
)
init_states = ComponentVector(
    SnowWater=0.0,
    LiquidWater=0.0,
    SoilWater=10.0,
    UpperZone=5.0,
    LowerZone=5.0,
)

model = DeepFlex.HBV(
    name="hbv",
    parameters=parameters,
    init_states=init_states
)

input = ComponentVector(Prcp=P, Pet=E, Temp=T)
output = DeepFlex.get_output(model, input=input, step=true)

# plot result
fig = Figure(size=(400, 300))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
x = range(1, 100, length=100)
lines!(ax, x, output[:Flow], color=:blue)
fig