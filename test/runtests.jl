using CSV
using DataFrames
using LumpedHydro
using LumpedHydro: step_func
using Test
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using Symbolics
using Lux
using LuxCore
using StableRNGs
using ComponentArrays
using DataInterpolations
using OrdinaryDiffEq
using Aqua

@testset "LumpedHydro.jl" begin
    include("run_flux.jl")
    include("utils/run_runtime_build.jl")
    include("run_bucket.jl")
    # Aqua.test_all(LumpedHydro; ambiguities = false)
end