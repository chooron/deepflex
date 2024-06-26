module LumpedHydro
## External packages
# common packages
using TOML
using Statistics
using Random
using ComponentArrays
using StructArrays
using NamedTupleTools
using Reexport
using StableRNGs
using DocStringExtensions
using BenchmarkTools

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
using Diffractor

# parameters Optimization
using Optimization
using OptimizationBBO
using OptimizationOptimisers

## package version
const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

## io type
const INPUT_TYPE = StructArray
const PAS_TYPE = ComponentVector
const FUNC_INPUT_TYPE = NamedTuple
const FUNC_PARAM_TYPE = NamedTuple

## Abstract Component Types
abstract type AbstractComponent end
abstract type AbstractSolver end
abstract type AbstractOptimizer end

abstract type AbstractEquation end

struct HydroEquation{input_names,output_names,param_names} <: AbstractEquation
    inputs::Vector{Num}
    outputs::Vector{Num}
    params::Vector{Num}

    function HydroEquation(
        input_names::Vector{Symbol},
        output_names::Vector{Symbol},
        param_names::Vector{Symbol},
    )
        inputs = vcat([@variables $var(t) = 0.0 for var in input_names]...)
        outputs = vcat([@variables $var(t) = 0.0 for var in output_names]...)
        params = vcat([@parameters $p = 0.0 [tunable = true] for p in param_names]...)
        return new{Tuple(input_names),Tuple(output_names),Tuple(param_names)}(inputs, outputs, params)
    end
end

abstract type AbstractFlux end
abstract type AbstractSimpleFlux <: AbstractFlux end
abstract type AbstractNeuralFlux <: AbstractFlux end
abstract type AbstractStateFlux <: AbstractFlux end
abstract type AbstractLagFlux <: AbstractFlux end

function (flux::AbstractFlux)(input::AbstractVector, params::AbstractVector)
    flux.inner_func(input, params)
end

function (flux::AbstractFlux)(input::AbstractMatrix, params::AbstractVector)
    flux.inner_func(input, params)
end

#* 负责某一平衡单元的计算
abstract type AbstractElement <: AbstractComponent end
#* 负责多个平衡联合单元的计算
abstract type AbstractUnit <: AbstractComponent end

# Sensealg type
const default_node_sensealg = BacksolveAdjoint(autojacvec=ZygoteVJP())
const default_ode_sensealg = ForwardDiffSensitivity()
# utils
include("utils/lossfunc.jl")
include("utils/name.jl")
include("utils/mtk.jl")
include("utils/optimize.jl")
include("utils/solver.jl")
include("utils/smoother.jl")
include("utils/graph.jl")
include("utils/unithydro.jl")
export step_func, ifelse_func

# framework build
include("flux.jl")
export SimpleFlux, StateFlux, LagFlux, NeuralFlux, StdMeanNormFlux, MinMaxNormFlux, TranparentFlux

include("element.jl")
export HydroElement, add_inputflux!, add_outputflux!, solve_prob

include("unit.jl")
export HydroUnit, update_unit!, add_elements!, remove_elements!

# Implement Flux
include("functions/evap.jl")
include("functions/flow.jl")
include("functions/infiltration.jl")
include("functions/melt.jl")
include("functions/others.jl")
include("functions/percolation.jl")
include("functions/pet.jl")
include("functions/rainfall.jl")
include("functions/recharge.jl")
include("functions/saturation.jl")
include("functions/snowfall.jl")
# Implements Models
# include("implements/PRNN.jl")
include("implements/cemaneige.jl")
include("implements/exphydro.jl")
include("implements/gr4j.jl")
include("implements/hbv_edu.jl")
include("implements/hymod.jl")
include("implements/simhyd.jl")
# include("implements/m50.jl")

# export abstract structs
export AbstractComponent, AbstractSolver, AbstractOptimizer,
    AbstractFlux, AbstractSimpleFlux, AbstractNeuralFlux, AbstractStateFlux, AbstractLagFlux,
    AbstractElement, AbstractUnit, AbstractNode

# export model
export ExpHydro, M50, GR4J, HyMOD, HBV

end