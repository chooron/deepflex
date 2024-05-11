module DeepFlex
## External packages
# common packages
using TOML
using Statistics
using Random
using ComponentArrays
using NamedTupleTools
using DataFrames
using Reexport
using StableRNGs

# run time stats
using BenchmarkTools

# ModelingToolkit building
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitNeuralNets
using ModelingToolkitStandardLibrary.Blocks
using Symbolics
using SymbolicUtils
using SciMLStructures: Tunable, replace, replace!

# graph compute
using Graphs
using MetaGraphs

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
abstract type AbstractParamInfo end
abstract type AbstractSmoother end
abstract type AbstractSolver end
abstract type AbstractOptimizer end

abstract type AbstractFlux end
abstract type AbstractSimpleFlux <: AbstractFlux end
abstract type AbstractNeuralFlux <: AbstractFlux end
abstract type AbstractStateFlux <: AbstractFlux end
abstract type AbstractLagFlux <: AbstractFlux end


abstract type AbstractReach end
#* 负责某一平衡单元的计算
abstract type AbstractElement <: AbstractComponent end
#* 负责多个平衡联合单元的计算
abstract type AbstractUnit <: AbstractComponent end
#* 负责单元的坡面汇流和洪水演进计算
abstract type AbstractNode <: AbstractComponent end
#* 负责基于节点的网络计算
abstract type AbstractNetwork <: AbstractComponent end
abstract type AbstractRiverNetwork <: AbstractNetwork end
abstract type AbstractGridNetwork <: AbstractNetwork end

# work for lux nn
Base.length(::Symbol) = 1

## Sensealg type
# const default_node_sensealg = BacksolveAdjoint(autojacvec=ZygoteVJP())
# const default_ode_sensealg = ForwardDiffSensitivity()
# utils
include("utils/graph.jl")
include("utils/lossfunc.jl")
include("utils/data.jl")
include("utils/mtk.jl")
include("utils/optimize.jl")
include("utils/solver.jl")
include("utils/smoother.jl")
include("utils/name.jl")
# framework build
include("flux.jl")
include("reach.jl")
include("element.jl")
include("node.jl")
include("network.jl")
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

# export
export SimpleFlux

end