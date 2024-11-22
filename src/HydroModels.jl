module HydroModels

## External packages
# common packages
using ComponentArrays
using ComponentArrays: indexmap, getval
using Dates
using DataFrames
using IterTools: ncycle
using LinearAlgebra
using NamedTupleTools
using ProgressMeter
using Reexport
using SparseArrays
using StableRNGs
using Statistics
using TOML

# runtime generated functions
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

# Symbolic building
using Symbolics
using Symbolics: tosymbol
using SymbolicUtils
using SymbolicUtils.Code
using ModelingToolkit: @variables, @parameters
using ModelingToolkit: t_nounits as t
using ModelingToolkit: isparameter
# graph compute
using Graphs

# data interpolataion
using DataInterpolations
using DataInterpolations: AbstractInterpolation

# solve ODEProblem
using SciMLBase
using OrdinaryDiffEq
using SciMLSensitivity

# deep learning
using Lux
using LuxCore
using NNlib

# parameters Optimization
using Optimization
using OptimizationBBO
using OptimizationOptimisers

## package version
const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])
const InputType = Union{AbstractArray,AbstractMatrix}

struct HydroEquation end
## Abstract Component Types
abstract type AbstractComponent end
abstract type AbstractHydroSolver end
abstract type AbstractHydroOptimizer end

abstract type AbstractFlux <: AbstractComponent end
abstract type AbstractElement <: AbstractComponent end
abstract type AbstractBucket <: AbstractElement end
abstract type AbstractRoute <: AbstractElement end
abstract type AbstractModel <: AbstractComponent end

abstract type AbstractHydroFlux <: AbstractFlux end
abstract type AbstractNeuralFlux <: AbstractHydroFlux end
abstract type AbstractStateFlux <: AbstractFlux end
abstract type AbstractUnitHydroFlux <: AbstractFlux end

abstract type AbstractDirectRoute <: AbstractRoute end
abstract type AbstractHydroRoute <: AbstractRoute end
abstract type AbstractRapidRoute <: AbstractRoute end

abstract type AbstractHydroWrapper <: AbstractComponent end
# utils
include("utils/uh.jl")
export UHFunction, UH_1_HALF, UH_2_FULL
include("utils/attr.jl")
include("utils/ca.jl")
include("utils/name.jl")
# include("utils/show.jl")
include("utils/build.jl")
include("utils/callback.jl")
include("utils/sort.jl")

include("optimizer.jl")
export BatchOptimizer, HydroOptimizer, GradOptimizer

include("solver.jl")
export ODESolver, DiscreteSolver, ManualSolver

# framework build
include("flux.jl")
export HydroFlux, StateFlux, NeuralFlux, UnitHydroFlux
include("bucket.jl")
export HydroBucket
include("route.jl")
export DirectRoute, GridRoute, VectorRoute, HydroRoute # , RapidRoute
include("model.jl")
export HydroModel
include("wrapper.jl")
export RecordComponentState, EstimateComponentParams, WeightSumComponentOutlet, ComputeComponentOutlet

# export abstract structs
export AbstractComponent, AbstractSolver, AbstractElement, AbstractUnit, AbstractHydroBucket, AbstractRoute, HydroEquation
export AbstractFlux, AbstractHydroFlux, AbstractNeuralFlux, AbstractStateFlux, AbstractRouteFlux

end # module HydroModels
