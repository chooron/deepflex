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
using LuxCUDA
using Zygote

# parameters Optimization
using Optimization
using OptimizationBBO
using Optimisers


# , DiffEqFlux, OrdinaryDiffEq, Optimization, OptimizationOptimJL,
#     OptimizationOptimisers, Random, Plots

## package version
const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])


## Abstract Component Types
abstract type Component end
abstract type AbstractUnit <: Component end
abstract type AbstractElement <: Component end
abstract type AbstractFunc <: Component end
abstract type AbstractNetwork <: Component end
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
# Optimization
include("framework/optimize.jl")
# Implements Models
include("implements/ExpHydro.jl")
# Implement Flux
include("functions/snowfall.jl")
include("functions/rainfall.jl")
include("functions/pet.jl")
include("functions/evap.jl")
include("functions/melt.jl")
include("functions/surfaceflow.jl")
include("functions/baseflow.jl")
# utils
include("utils/smooth_func.jl")
include("utils/loss_func.jl")
include("utils/copy.jl")
include("utils/graph_utils.jl")
end