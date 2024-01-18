module DeepFlex
## External packages
# common packages
using Parameters
using TOML

# graph compute
using Graphs
using MetaGraphs

# data interpolataion
using Interpolations

# solve ODEProblem
using DifferentialEquations

## package version
const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

## Abstract Element Types
abstract type BaseElement end
abstract type ParameterizedElement <: BaseElement end
abstract type StateElement <: BaseElement end
abstract type StateParameterizedElement <: BaseElement end
abstract type ODEsElement <: StateParameterizedElement end
abstract type DiscElement <: StateParameterizedElement end

## Abstract Unit Type
abstract type Component end
abstract type Unit end

export ODEsElement

# Element Methods
include("framework/element.jl")
# Implements Element
include("elements/exphydro.jl")
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
end