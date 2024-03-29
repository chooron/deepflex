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
using Symbolics

# graph compute
using Graphs
using MetaGraphs

# data interpolataion
using DataInterpolations
# solve ODEProblem
using OrdinaryDiffEq
using DiffEqFlux

# deep learning
using Lux
using Zygote

# parameters Optimization
using Optimization
using OptimizationBBO
using OptimizationOptimisers

## package version
const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

## Abstract Component Types
abstract type AbstractComponent end

abstract type AbstractParamInfo end
abstract type AbstractSmoother end
abstract type AbstractSolver end
abstract type AbstractOptimizer end
abstract type AbstractFlux end
abstract type AbstractElement end

abstract type AbstractUnit <: AbstractComponent end
abstract type AbstractNode <: AbstractComponent end
abstract type AbstractNetwork <: AbstractComponent end

## Sensealg type
const default_node_sensealg = BacksolveAdjoint(autojacvec=ZygoteVJP())
const default_ode_sensealg = ForwardDiffSensitivity()

## for ModelingToolkit
@variables t
const D = Differential(t)

# framework build
include("framework/solver.jl")
include("framework/fluxes.jl")
include("framework/element.jl")
include("framework/unit.jl")
include("framework/node.jl")
include("framework/network.jl")
# framework methods
include("framework/optimize.jl")
# Implement Flux
include("functions/baseflow.jl")
include("functions/evap.jl")
include("functions/flow.jl")
include("functions/infiltration.jl")
include("functions/melt.jl")
include("functions/miscellaneous.jl")
include("functions/percolation.jl")
include("functions/pet.jl")
include("functions/rainfall.jl")
include("functions/recharge.jl")
include("functions/saturation.jl")
include("functions/smoother.jl")
include("functions/snowfall.jl")
include("functions/surfaceflow.jl")
include("functions/unithydro.jl")
# Implement Element
include("elements/slope.jl")
include("elements/soil.jl")
include("elements/surface.jl")
# Implements Models
include("implements/ExpHydro.jl")
include("implements/GR4J.jl")
include("implements/HyMOD.jl")
include("implements/HBV.jl")
include("implements/XAJ.jl")
include("implements/HydroNODE.jl")
# utils
include("utils/loss_func.jl")
include("utils/graph_utils.jl")
# Optimization
include("framework/optimize.jl")
end