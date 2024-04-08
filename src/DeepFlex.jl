module DeepFlex
## External packages
# common packages
using TOML
using Statistics
using Random
using ComponentArrays
using NamedTupleTools
using DataFrames

# run time stats
using BenchmarkTools

# ModelingToolkit building
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Symbolics

# graph compute
using Graphs
using MetaGraphs

# data interpolataion
using DataInterpolations
# solve ODEProblem
using OrdinaryDiffEq
# using DiffEqFlux

# deep learning
using Lux
using Zygote

# parameters Optimization
using Optimization
using OptimizationBBO
using OptimizationOptimisers

## package version
# const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

## Abstract Component Types
abstract type AbstractComponent end

abstract type AbstractParamInfo end
abstract type AbstractSmoother end
abstract type AbstractSolver end
abstract type AbstractOptimizer end

abstract type AbstractFlux end
abstract type AbstractNNFlux <: AbstractFlux end

#* 负责某一平衡单元的计算
abstract type AbstractElement <: AbstractComponent end
#* 负责多个平衡联合单元的计算
abstract type AbstractUnit <: AbstractComponent end
#* 负责单元的坡面汇流和洪水演进计算
abstract type AbstractNode <: AbstractComponent end
#* 负责单元的河网汇流计算
abstract type AbstractNetwork <: AbstractComponent end

## Sensealg type
# const default_node_sensealg = BacksolveAdjoint(autojacvec=ZygoteVJP())
# const default_ode_sensealg = ForwardDiffSensitivity()
# utils
include.(filter(contains(r".jl$"), readdir("utils"; join=true)))
# framework build
include("flux.jl")
include("element.jl")
include("unit.jl")
include("node.jl")
include("network.jl")
# Implement Flux
include.(filter(contains(r".jl$"), readdir("functions"; join=true)))
# Implement Element
# Implements Models
include.(filter(contains(r".jl$"), readdir("implements"; join=true)))
end