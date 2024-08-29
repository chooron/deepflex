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
    # include("run_flux.jl")
    include("run_bucket.jl")
    # Aqua.test_all(LumpedHydro; ambiguities = false)
end