using CSV
using DataFrames
using Lux
using Test
using ModelingToolkit
using Symbolics
using LuxCore
using ComponentArrays
using DataInterpolations
using OrdinaryDiffEq
using Statistics
using Graphs
# using HydroModels
include("../src/HydroModels.jl")

@variables q_in q_out s_river q_gen
@parameters lag

exprs = [(s_river + q_in) / (1 + lag), q_in - (s_river + q_in) / (1 + lag)]
rflux = HydroModels.StateRouteFlux(q_in, q_out, [lag], s_river, exprs=exprs)
rflux.func(10.0, 0.0, 0.2, 1.0)