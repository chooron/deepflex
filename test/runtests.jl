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
using Graphs
using HydroModels
using Test

@testset "HydroModels.jl" begin
    include("run_flux.jl")
    include("run_bucket.jl")
    include("run_uh.jl")
    include("run_route.jl")
    include("run_lumped_model.jl")
    include("run_spatial_model.jl")
end