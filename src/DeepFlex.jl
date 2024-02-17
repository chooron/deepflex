module DeepFlex
## External packages
# common packages
using TOML
using Statistics
using Random
using ComponentArrays

# graph compute
using Graphs
using MetaGraphs

# data interpolataion
using Interpolations

# solve ODEProblem
using OrdinaryDiffEq
using DiffEqFlux

# deep learning
using Lux
using Zygote

# parameters Optimization
using Optimization
using OptimizationBBO
# using Optimisers
using OptimizationOptimisers


# , DiffEqFlux, OrdinaryDiffEq, Optimization, OptimizationOptimJL,
#     OptimizationOptimisers, Random, Plots

## package version
const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])
const ode_solve = "ode_solve"
const dct_solve = "dct_solve"

## Abstract Component Types
abstract type AbstractComponent end
abstract type AbstractUnit <: AbstractComponent end
abstract type AbstractElement <: AbstractComponent end
abstract type AbstractFunc <: AbstractComponent end
abstract type SimpleFunc <: AbstractFunc end
abstract type LuxNNFunc <: AbstractFunc end
abstract type AbstractNetwork <: AbstractComponent end


## Sensealg type
const default_node_sensealg = BacksolveAdjoint(autojacvec=ZygoteVJP())
const default_ode_sensealg = ForwardDiffSensitivity()

## Abstract Element Types
# abstract type LagElement <: AbstractElement end
# abstract type ODEElement <: AbstractElement end
# abstract type DiscElement <: AbstractElement end
# abstract type LuxElement <: StateParameterizedElement end

# framework Methods
include("framework/paraminfo.jl")
include("framework/element.jl")
include("framework/unit.jl")
include("framework/node.jl")
include("framework/network.jl")
# Implement Element
include("elements/routingstore.jl")
include("elements/snowwater.jl")
include("elements/soilWater.jl")
# Implement Flux
include("functions/baseflow.jl")
include("functions/evap.jl")
include("functions/flow.jl")
include("functions/melt.jl")
include("functions/nn.jl")
include("functions/percolation.jl")
include("functions/pet.jl")
include("functions/rainfall.jl")
include("functions/recharge.jl")
include("functions/routing.jl")
include("functions/saturation.jl")
include("functions/snowfall.jl")
include("functions/splitter.jl")
include("functions/surfaceflow.jl")
include("functions/transparent.jl")
# Implements Models
include("implements/ExpHydro.jl")
include("implements/GR4J.jl")
include("implements/HydroNODE.jl")
# utils
include("utils/smooth_func.jl")
include("utils/loss_func.jl")
include("utils/graph_utils.jl")
# Optimization
include("framework/optimize.jl")
end