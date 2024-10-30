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
using SymbolicUtils
using SymbolicUtils.Code
using ModelingToolkit: @variables, @parameters
using ModelingToolkit: t_nounits as t
# graph compute
using Graphs

# data interpolataion
using DataInterpolations
using DataInterpolations: AbstractInterpolation

# integral
using Integrals

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
# HydroErrors
using HydroErrors

# HydroEquations
# using HydroEquations

## package version
const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

## Abstract Component Types
abstract type AbstractComponent end
## todo this is used for multiple propose running
abstract type AbstractRunner end
abstract type AbstractSolver end

#* 负责某一平衡单元的计算
abstract type AbstractFlux <: AbstractComponent end
abstract type AbstractElement <: AbstractComponent end
abstract type AbstractBucket <: AbstractElement end
abstract type AbstractRoute <: AbstractElement end
#* 负责多个平衡联合单元的计算
abstract type AbstractModel <: AbstractComponent end

#* 参数和状态的估计
abstract type AbstractEstimator <: AbstractComponent end

#* Flux的多种变体
abstract type AbstractSimpleFlux <: AbstractFlux end
abstract type AbstractNeuralFlux <: AbstractFlux end
abstract type AbstractStateFlux <: AbstractFlux end
abstract type AbstractRouteFlux <: AbstractFlux end
abstract type AbstractStateRouteFlux <: AbstractRouteFlux end
abstract type AbstractDynamicRouteFlux <: AbstractRouteFlux end
abstract type AbstractUnitHydroFlux <: AbstractFlux end
abstract type AbstractTimeVaryingFlux <: AbstractFlux end

#* route的两种变体
abstract type AbstractSumRoute <: AbstractRoute end
abstract type AbstractHydroRoute <: AbstractRoute end
abstract type AbstractRapidRoute <: AbstractRoute end

"""
Metadata about the component, including:
- name: Symbol representing the component's name
- inputs: Vector of input variable names
- outputs: Vector of output variable names
- params: Vector of parameter names
- states: Vector of state variable names
- nns: Vector of neural network names
"""
struct HydroMeta
    name::Symbol
    inputs::Vector{Symbol}
    outputs::Vector{Symbol}
    params::Vector{Symbol}
    states::Vector{Symbol}
    nns::Vector{Symbol}

    function HydroMeta(;
        name::Symbol, inputs::Vector{Symbol}=Symbol[], outputs::Vector{Symbol}=Symbol[],
        params::Vector{Symbol}=Symbol[], states::Vector{Symbol}=Symbol[], nns::Vector{Symbol}=Symbol[],
    )
        return new(name, inputs, outputs, params, states, nns)
    end
end

# utils
include("utils/attr.jl")
include("utils/ca.jl")
include("utils/name.jl")
include("utils/show.jl")
include("utils/build.jl")
include("utils/callback.jl")
include("utils/sort.jl")

include("optimize.jl")
export param_grad_optim, param_box_optim, nn_param_optim, param_batch_optim

include("solver.jl")
export ODESolver, DiscreteSolver, ManualSolver

# framework build
include("flux.jl")
export SimpleFlux, StateFlux, NeuralFlux, HydroFunction
include("uhflux.jl")
export UHFunction, UH_1_HALF, UH_2_FULL, UnitHydroFlux
include("rflux.jl")
export RouteFlux, RiverRouteFlux
include("bucket.jl")
export HydroBucket
include("route.jl")
export WeightSumRoute, GridRoute, VectorRoute, RapidRoute
include("model.jl")
export HydroModel
include("estimator.jl")
export HydroEstimator

# export abstract structs
export AbstractComponent, AbstractSolver, AbstractElement, AbstractUnit, AbstractHydroBucket, AbstractRoute
export AbstractFlux, AbstractSimpleFlux, AbstractNeuralFlux, AbstractStateFlux, AbstractRouteFlux

# export model
export ExpHydro, M50, GR4J, HyMOD, HBV_EDU

end # module HydroModels
