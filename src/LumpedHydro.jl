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
using ModelingToolkitNeuralNets
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
        if vk in vks2
            _p = vcat(_p, ca2[vk])
        else
            _p = vcat(_p, ca[vk])
        end
    end
    ComponentArray(_p, ax)
end

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

    function HydroEquation(
        inputs::Vector{Num},
        outputs::Vector{Num},
        params::Vector{Num},
    )
        input_names = Symbolics.tosymbol.(inputs, escape=false)
        output_names = Symbolics.tosymbol.(outputs, escape=false)
        param_names = Symbolics.tosymbol.(params, escape=false)
        return new{Tuple(input_names),Tuple(output_names),Tuple(param_names)}(inputs, outputs, params)
    end
end

abstract type AbstractFlux <: AbstractComponent end
abstract type AbstractSimpleFlux <: AbstractFlux end
abstract type AbstractNeuralFlux <: AbstractSimpleFlux end
abstract type AbstractStateFlux <: AbstractFlux end
abstract type AbstractLagFlux <: AbstractFlux end

function (flux::AbstractFlux)(input::AbstractArray, params::AbstractArray)
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
include("utils/show.jl")
include("utils/ode.jl")
include("utils/solver.jl")
include("utils/special_fluxes.jl")
include("utils/smoother.jl")
include("utils/graph.jl")
include("utils/unithydro.jl")
export step_func, ifelse_func

# framework build
include("flux.jl")
export SimpleFlux, StateFlux, LagFlux, NeuralFlux

include("element.jl")
export HydroElement, solve_prob # , add_inputflux!, add_outputflux!, 

include("unit.jl")
export HydroUnit #, update_unit!, add_elements!, remove_elements!

include("optimize.jl")
export param_grad_optim, param_box_optim

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
# include("implements/PRNN.jl")
include("implements/cemaneige.jl")
include("implements/exphydro.jl")
include("implements/gr4j.jl")
include("implements/hbv_edu.jl")
include("implements/hymod.jl")
include("implements/simhyd.jl")
include("implements/m50.jl")

# export abstract structs
export AbstractComponent, AbstractSolver, AbstractOptimizer,
    AbstractFlux, AbstractSimpleFlux, AbstractNeuralFlux, AbstractStateFlux, AbstractLagFlux,
    AbstractElement, AbstractUnit, AbstractNode

# export model
export ExpHydro, M50, GR4J, HyMOD, HBV

# export special flux
export StdMeanNormFlux, MinMaxNormFlux, TranparentFlux

end