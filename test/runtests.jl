using Aqua
using CSV
using DataFrames
using Lux
using Test
using ModelingToolkit
using Symbolics
using LuxCore
using StableRNGs
using ComponentArrays
using DataInterpolations
using OrdinaryDiffEq
using Statistics
include("../src/HydroModels.jl")

@testset "HydroModels.jl" begin
    include("run_bucket.jl")
    # include("run_route.jl")
    # Aqua.test_all(LumpedHydro; ambiguities = false)
end