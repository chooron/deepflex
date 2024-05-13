using ComponentArrays
using ModelingToolkit
using DifferentialEquations

include("../../src/LumpedHydro.jl")

tmp_input = Float64[1, 2, 3, 4, 5, 3, 2, 1, 4]

rf = LumpedHydro.LagFlux(:q, :ql, lag_func=LumpedHydro.uh_1_half, param_names=:x4)

rf((q=tmp_input,), (x4=3.5,))