module DeepFlex
## External packages
# common packages
using TOML
using Statistics

# graph compute
using Graphs
using MetaGraphs

# data interpolataion
using Interpolations

# solve ODEProblem
using DifferentialEquations

# parameters Optimizationusing 
using Optimization
using OptimizationBBO

## package version
const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

## Abstract Element Types
abstract type BaseElement end
abstract type ParameterizedElement <: BaseElement end
abstract type StateElement <: BaseElement end
abstract type StateParameterizedElement <: BaseElement end
abstract type LagElement <: StateParameterizedElement end
abstract type ODEsElement <: StateParameterizedElement end
abstract type DiscElement <: StateParameterizedElement end
# abstract type LuxElement <: StateParameterizedElement end

## Abstract Component Types
abstract type Component end
abstract type Unit <: Component end

export ODEsElement


# framework Methods
include("framework/paraminfo.jl")
include("framework/element.jl")
include("framework/unit.jl")
include("framework/node.jl")
include("framework/network.jl")
# Optimization
include("framework/optimize.jl")
# Implements Models
include("models/exphydro.jl")
# Implement Flux
include("fluxes/snowfall.jl")
include("fluxes/rainfall.jl")
include("fluxes/pet.jl")
include("fluxes/evap.jl")
include("fluxes/melt.jl")
include("fluxes/surfaceflow.jl")
include("fluxes/baseflow.jl")
# utils
include("utils/smooth_func.jl")
include("utils/loss_func.jl")
include("utils/copy.jl")
include("utils/graph_utils.jl")
end