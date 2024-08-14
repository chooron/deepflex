module LumpedHydro
## External packages
# common packages
using TOML
using Dates
using Reexport
using StableRNGs
using Statistics
using ComponentArrays
using ComponentArrays: indexmap, getval
using NamedTupleTools
using DocStringExtensions
using IterTools: ncycle
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

# ModelingToolkit building
using ModelingToolkit
using ModelingToolkit: D_nounits as D
using Symbolics
using SymbolicUtils
using SymbolicUtils.Code

# graph compute
using Graphs

# data interpolataion
using DataInterpolations
# using Interpolations

# solve ODEProblem
using SciMLBase
using OrdinaryDiffEq
using DiffEqFlux

# deep learning
using Lux
using Zygote

# parameters Optimization
using Optimization
using OptimizationBBO
using OptimizationOptimisers

# HydroErr
using HydroErr

# HydroEquations
# using HydroEquations

## package version
const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

## Abstract Component Types
abstract type AbstractComponent end
abstract type AbstractSolver end

abstract type AbstractFlux <: AbstractComponent end
abstract type AbstractSimpleFlux <: AbstractFlux end
abstract type AbstractNeuralFlux <: AbstractSimpleFlux end
abstract type AbstractStateFlux <: AbstractFlux end
abstract type AbstractLagFlux <: AbstractFlux end

#* 负责某一平衡单元的计算
abstract type AbstractElement <: AbstractComponent end
abstract type AbstractHydroElement <: AbstractElement end
abstract type AbstractLagElement <: AbstractElement end
abstract type AbstractRoute <: AbstractElement end
#* 负责多个平衡联合单元的计算
abstract type AbstractUnit <: AbstractComponent end

# Sensealg type
const default_node_sensealg = BacksolveAdjoint(autojacvec=ZygoteVJP())
const default_ode_sensealg = ForwardDiffSensitivity()
# utils
include("utils/attr.jl")
include("utils/ca.jl")
include("utils/name.jl")
include("utils/show.jl")
include("utils/graph.jl")
include("utils/smooth.jl")
include("utils/unithydro.jl")

# framework build
include("flux.jl")
export SimpleFlux, StateFlux, LagFlux, NeuralFlux
# special flux
include("utils/normalize.jl")
export StdMeanNormFlux, MinMaxNormFlux, TranparentFlux

include("element.jl")
export HydroElement, LagElement, solve_prob # , add_inputflux!, add_outputflux!, 

include("unit.jl")
export HydroUnit #, update_unit!, add_elements!, remove_elements!

include("optimize.jl")
export param_grad_optim, param_box_optim, nn_param_optim

include("solver.jl")
export ODESolver, DiscreteSolver, ManualSolver

# Implements Models
include("implements/cemaneige.jl")
include("implements/exphydro.jl")
include("implements/gr4j.jl")
include("implements/hbv_edu.jl")
include("implements/hbv_maxbas.jl")
include("implements/hbv_nn.jl")
include("implements/hymod.jl")
include("implements/simhyd.jl")
include("implements/m50.jl")

# export abstract structs
export AbstractComponent, AbstractSolver, AbstractElement, AbstractUnit
export AbstractFlux, AbstractSimpleFlux, AbstractNeuralFlux, AbstractStateFlux, AbstractLagFlux

# export model
export ExpHydro, M50, GR4J, HyMOD, HBV_EDU
end