module DeepFlex
# External packages
using CSV
using Parameters
using DifferentialEquations
using TOML


const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])
# Abstract Element Types
abstract type BaseElement end
abstract type ParameterizedElement <: BaseElement end
abstract type StateElement <: BaseElement end
abstract type StateParameterizedElement <: BaseElement end
abstract type ODEsElement <: StateParameterizedElement end
abstract type DiscElement <: StateParameterizedElement end

# Abstract Unit Type
abstract type Component end

# Abstract Flux Type
abstract type Flux end


# Element Methods
include("framework/element.jl")


# Implements Element
include("elements/exphydro.jl")

# Implement Flux
include("fluxes/snowfall.jl")

end