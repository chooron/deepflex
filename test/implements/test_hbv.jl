# import lib
using CSV
using Random
using DataFrames
using CairoMakie
using ComponentArrays
using CairoMakie: Axis

# test gr4j model
include("../../src/LumpedHydro.jl")

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

unit_params = ComponentVector(
    tt=tt, tti=tti, cfr=cfr, cfmax=cfmax, ttm=ttm, whc=whc,
    cflux=cflux, fc=fc, lp=lp, k0=k0, k1=k1, α=α, β=β, c=c, maxbas=maxbas
)
init_states = ComponentVector(
    snowwater=0.0,
    liquidwater=0.0,
    soilwater=10.0,
    upperzone=5.0,
    lowerzone=5.0,
)

ps = ComponentVector(hbv=(params=unit_params, initstates=init_states, weight=1.0))

model = LumpedHydro.HBV.Node(name=:hbv, mtk=true, step=true);
LumpedHydro.get_param_names(model.units[1].surface[1])
input = (hbv=(prcp=P, pet=E, temp=T, time=1:length(T)),)
output = model(input, ps)

# plot result
fig = Figure(size=(400, 300))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
x = range(1, 100, length=100)
lines!(ax, x, output[:flow], color=:blue)
fig