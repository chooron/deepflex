module LumpedHydro
## External packages
# common packages
using TOML
using Reexport
using StableRNGs
using Statistics
using ComponentArrays
using NamedTupleTools
using DocStringExtensions
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

# ModelingToolkit building
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Symbolics
using SymbolicUtils
using SymbolicUtils.Code

# graph compute
using Graphs

# data interpolataion
using DataInterpolations

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

## package version
const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

## Abstract Component Types
abstract type AbstractComponent end
abstract type AbstractSolver end
abstract type AbstractEquation end

# merge ComponentArray:https://github.com/jonniedie/ComponentArrays.jl/issues/186
function merge_ca(ca::ComponentArray{T1}, ca2::ComponentArray{T2}) where {T1,T2}
    ax = getaxes(ca)
    ax2 = getaxes(ca2)
    vks = valkeys(ax[1])
    vks2 = valkeys(ax2[1])
    _p = Vector{T2}()
    for vk in vks
        if length(getaxes(ca[vk])) > 0
            _p = vcat(_p, collect(merge_ca(ca[vk], vk in vks2 ? getproperty(ca2, vk) : ComponentVector())))
        else
            if vk in vks2
                _p = vcat(_p, ca2[vk])
            else
                _p = vcat(_p, ca[vk])
            end
        end
    end
    ComponentArray(_p, ax)
end

abstract type AbstractFlux <: AbstractComponent end
abstract type AbstractSimpleFlux <: AbstractFlux end
abstract type AbstractNeuralFlux <: AbstractSimpleFlux end
abstract type AbstractStateFlux <: AbstractFlux end
abstract type AbstractLagFlux <: AbstractFlux end

(flux::AbstractFlux)(input::AbstractArray, params::AbstractArray) = flux.inner_func(input, params)

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
include("utils/show.jl")
include("utils/smoother.jl")
include("utils/graph.jl")
include("utils/unithydro.jl")
export step_func, ifelse_func

# framework build
include("equation.jl")

include("flux.jl")
export SimpleFlux, StateFlux, LagFlux, NeuralFlux
# special flux
include("utils/normalize.jl")
export StdMeanNormFlux, MinMaxNormFlux, TranparentFlux

include("element.jl")
export HydroElement, solve_prob # , add_inputflux!, add_outputflux!, 

include("unit.jl")
export HydroUnit #, update_unit!, add_elements!, remove_elements!

include("optimize.jl")
export param_grad_optim, param_box_optim, nn_param_optim

include("solver.jl")
export ODESolver, DiscreteSolver, ManualSolver

# Implement Flux
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
# Implements Models
include("implements/cemaneige.jl")
include("implements/exphydro.jl")
include("implements/gr4j.jl")
include("implements/hbv_edu.jl")
include("implements/hymod.jl")
include("implements/simhyd.jl")
include("implements/m50.jl")

# export abstract structs
export AbstractComponent, AbstractSolver, AbstractElement, AbstractUnit
export AbstractFlux, AbstractSimpleFlux, AbstractNeuralFlux, AbstractStateFlux, AbstractLagFlux

# export model
export ExpHydro, M50, GR4J, HyMOD, HBV_EDU
end