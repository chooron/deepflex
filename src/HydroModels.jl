module HydroModels

## External packages
# common packages
using TOML
using Dates
using Reexport
using StableRNGs
using Statistics
using SparseArrays
using ComponentArrays
using ComponentArrays: indexmap, getval
using NamedTupleTools
using DocStringExtensions
using IterTools: ncycle
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

# solve ODEProblem
using SciMLBase
using OrdinaryDiffEq
using SciMLSensitivity

# deep learning
using Lux
using LuxCore
using NNlib
using Zygote

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
abstract type AbstractUnitHydroFlux <: AbstractFlux end
abstract type AbstractTimeVaryingFlux <: AbstractFlux end

#* route的两种变体
abstract type AbstractSumRoute <: AbstractRoute end
abstract type AbstractGridRoute <: AbstractRoute end
abstract type AbstractVectorRoute <: AbstractRoute end

# utils
include("utils/attr.jl")
include("utils/ca.jl")
include("utils/name.jl")
include("utils/show.jl")
include("utils/graph.jl")
include("utils/smooth.jl")
include("utils/build.jl")
# some unit hydro function
include("utils/unithydro.jl")

# framework build
include("flux.jl")
export SimpleFlux, StateFlux, NeuralFlux, RouteFlux, UnitHydroFlux

include("bucket.jl")
export HydroBucket # , add_inputflux!, add_outputflux!, 

include("route.jl")
export WeightSumRoute, GridRoute, VectorRoute

include("model.jl")
export HydroModel #, update_unit!, add_elements!, remove_elements!

include("estimator.jl")
export HydroEstimator

include("optimize.jl")
export param_grad_optim, param_box_optim, nn_param_optim

include("solver.jl")
export ODESolver, DiscreteSolver, ManualSolver

# some route function and special flux
include("fluxes/cascade.jl")
export CascadeRouteFlux
include("fluxes/discharge.jl")
export DischargeRouteFlux
include("fluxes/muskingum.jl")
export MuskingumRouteFlux
include("fluxes/normalize.jl")
export StdMeanNormFlux, MinMaxNormFlux
include("fluxes/rename.jl")
export RenameFlux

# Implements Models
include("implements/cemaneige.jl")
include("implements/exphydro.jl")
include("implements/gr4j.jl")
include("implements/hymod.jl")
include("implements/simhyd.jl")
include("implements/m50.jl")

# export abstract structs
export AbstractComponent, AbstractSolver, AbstractElement, AbstractUnit, AbstractHydroBucket, AbstractRoute
export AbstractFlux, AbstractSimpleFlux, AbstractNeuralFlux, AbstractStateFlux, AbstractRouteFlux

# export model
export ExpHydro, M50, GR4J, HyMOD, HBV_EDU

end # module HydroModels
