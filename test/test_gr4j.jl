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
x1, x2, x3, x4 = 50.0, 0.1, 20.0, 3.5
parameters = ComponentVector(x1=x1, x2=x2, x3=x3, x4=x4, ω=3.5, γ=5.0)
init_states = ComponentVector(SoilWater=0.0, RoutingStore=10.0)
model = DeepFlex.GR4J(
    name="gr4j",
    parameters=parameters,
    init_states=init_states
)
input = ComponentVector(Prcp=P, Pet=E)

output = DeepFlex.get_output(model, input=input, step=true)