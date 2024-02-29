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
P[1:10] = rand(0.0:0.01:10.0, 10)
P[26:30] = rand(0.0:0.01:20.0, 5)
P[41:60] = rand(0.0:0.01:5.0, 20)
P[81:83] = rand(30.0:0.01:50.0, 3)

# build model
aim, b, a, stot, fwm, flm, c, ex, ki, kg, ci, cg = 0.6, 0.1, 5, 500, 0.5, 0.5, 0.5, 5, 0.5, 0.5, 0.5, 0.5
parameters = ComponentVector(Aim=aim, Wmax=fwm * stot, Smax=(1 - fwm) * stot,
    b=b, a=a, c=c, LM=flm * fwm * stot,
    stot=stot, fwm=fwm, flm=flm,
    ex=ex, ki=ki, kg=kg, ci=ci, cg=cg)

init_states = ComponentVector(
    TensionWater=10.0, FreeWater=10.0,
    InterRouting=5.0, BaseRouting=5.0
)

model = DeepFlex.XAJ(
    name="xaj",
    parameters=parameters,
    init_states=init_states
)

input = ComponentVector(Prcp=P, Pet=E)
output = DeepFlex.get_output(model, input=input, step=true)

# plot result
fig = Figure(size=(400, 300))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
x = range(1, 100, length=100)
lines!(ax, x, output[:Flow], color=:blue)
fig