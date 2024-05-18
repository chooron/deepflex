module LumpedHydro
## External packages
# common packages
using TOML
using Statistics
using Random
using ComponentArrays
using NamedTupleTools
using Reexport
using StableRNGs
using DocStringExtensions

# Multitreading and parallel computing
using Base.Threads

# ModelingToolkit building
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitNeuralNets
using ModelingToolkitStandardLibrary.Blocks
using Symbolics
using SymbolicUtils
using SymbolicIndexingInterface
using SciMLStructures: Tunable, replace

# graph compute
using Graphs

# data interpolataion
using DataInterpolations

# solve ODEProblem
using SciMLBase
using OrdinaryDiffEq
using DiffEqFlux

# solver NonlinearProblem
using NonlinearSolve

# deep learning
using Lux
using Zygote

# parameters Optimization
using Optimization
using OptimizationBBO
using OptimizationOptimisers

## package version
const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

## Abstract Component Types
abstract type AbstractComponent end
abstract type AbstractSolver end
abstract type AbstractOptimizer end

abstract type AbstractFlux end
abstract type AbstractSimpleFlux <: AbstractFlux end
abstract type AbstractNeuralFlux <: AbstractFlux end
abstract type AbstractStateFlux <: AbstractFlux end
abstract type AbstractLagFlux <: AbstractFlux end

#* 负责某一平衡单元的计算
abstract type AbstractElement <: AbstractComponent end
#* 负责多个平衡联合单元的计算
abstract type AbstractUnit <: AbstractComponent end
#* 负责单元的内部汇流计算
abstract type AbstractNode <: AbstractComponent end

# work for lux nn
Base.length(::Symbol) = 1

## Sensealg type
# const default_node_sensealg = BacksolveAdjoint(autojacvec=ZygoteVJP())
# const default_ode_sensealg = ForwardDiffSensitivity()
# utils
include("utils/lossfunc.jl")
include("utils/data.jl")
include("utils/mtk.jl")
include("utils/optimize.jl")
include("utils/solver.jl")
include("utils/smoother.jl")
include("utils/name.jl")
include("utils/graph.jl")

# framework build
include("flux.jl")
export SimpleFlux, StateFlux, LagFlux, NeuralFlux, StdMeanNormFlux, MinMaxNormFlux, TranparentFlux

include("element.jl")
export HydroElement, add_inputflux!, add_outputflux!, solve_prob

include("unit.jl")
export HydroUnit, update_unit!, add_elements!, remove_elements!

include("node.jl")
export HydroNode

# Implement Flux
include("functions/baseflow.jl")
include("functions/evap.jl")
include("functions/flow.jl")
include("functions/infiltration.jl")
include("functions/melt.jl")
include("functions/percolation.jl")
include("functions/pet.jl")
include("functions/rainfall.jl")
include("functions/recharge.jl")
include("functions/saturation.jl")
include("functions/snowfall.jl")
include("functions/surfaceflow.jl")
include("functions/unithydro.jl")
# Implements Models
# include("implements/PRNN.jl")
include("implements/exphydro.jl")
include("implements/m50.jl")
include("implements/gr4j.jl")
include("implements/hymod.jl")
include("implements/hbv.jl")

# export abstract structs
export AbstractComponent, AbstractSolver, AbstractOptimizer,
    AbstractFlux, AbstractSimpleFlux, AbstractNeuralFlux, AbstractStateFlux, AbstractLagFlux,
    AbstractElement, AbstractUnit, AbstractNode

# export model
export ExpHydro, M50, GR4J, HyMOD, HBV

end