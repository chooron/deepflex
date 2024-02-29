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
abstract type AbstractParamInfo <: AbstractComponent end
abstract type AbstractFlux <: AbstractComponent end
abstract type AbstractUnit <: AbstractComponent end
abstract type AbstractElement <: AbstractComponent end
abstract type AbstractNetwork <: AbstractComponent end


## Sensealg type
const default_node_sensealg = BacksolveAdjoint(autojacvec=ZygoteVJP())
const default_ode_sensealg = ForwardDiffSensitivity()

# framework Methods
include("framework/paraminfo.jl")
include("framework/smooth.jl")
include("framework/fluxes.jl")
include("framework/element.jl")
include("framework/unit.jl")
include("framework/node.jl")
include("framework/network.jl")
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
include("functions/snowfall.jl")
include("functions/surfaceflow.jl")
include("functions/unithydro.jl")
# Implement Element
include("elements/lag.jl")
include("elements/routing.jl")
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
include("utils/smooth_func.jl")
include("utils/loss_func.jl")
include("utils/graph_utils.jl")
# Optimization
include("framework/optimize.jl")
end