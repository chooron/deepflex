module Element

# 定义一个虚拟的Element类
abstract type BaseElement end
abstract type ParameterizedElement <: BaseElement end
abstract type StateElement <: BaseElement end
abstract type StateParameterizedElement <: BaseElement end
abstract type ODEsElement <: StateParameterizedElement end
abstract type DiscElement <: StateParameterizedElement end


function get_id(ele::BaseElement)::String
    return ele.id
end

function get_parameters(ele::Union{ParameterizedElement,StateParameterizedElement}; names::Vector{String}=nothing)::Dict{String, Any}
    if isnothing(names)
        return ele.parameters
    else
        return Dict(name => ele.parameters[name] for name in names)
    end
end

function get_states(ele::Union{StateElement,StateParameterizedElement}; names::Vector{String}=nothing)::Dict{String, Vector{Number}}
    if isnothing(names)
        return ele.states
    else
        return Dict(name => ele.states[name] for name in names)
    end
end

function solve_prob(ele::ODEsElement, solved::bool, solver::Any)
    if !solved
        """
        Todo: use ode solver
        """
        state_arr = solver.solve(ele.flux_funs, ele.init_states, )
    end
end
end
 
# mutable struct BaseElement <: BaseElementType
#     id::String
# end

# function BaseElement(id::String)
#     BaseElement(id)
# end

# mutable struct ParameterizedElement <: ParameterizedElementType
#     id::String
#     parameters::Dict{String,Any}
# end

# function ParameterizedElement(id::String, parameters::Dict{String,Any})
#     ParameterizedElement(id, parameters)
# end


# mutable struct StateElement <: StateElementType
#     id::String
#     init_states::Dict{String,Number}

#     # 中间计算状态
#     states::Dict{String,Vector{Number}}
# end

# function StateElement(id::String, states::Dict{String,Number})
#     StateElement(id, states, nothing)
# end


# mutable struct StateParameterizedElement <: StateParameterizedElementType
#     id::String
#     parameters::Dict{String,Any}
#     init_states::Dict{String,Number}

#     # 中间计算状态
#     states::Dict{String,Vector{Number}}
# end

# function StateParameterizedElement(id::String, parameters::Dict{String,Any}, states::Dict{String,Number})
#     StateParameterizedElement(id, parameters, states, nothing)
# end


# function get_parameters(ele::Union{ParameterizedElement,StateParameterizedElement}; names::Vector{String}=nothing)
#     if isnothing(names)
#         return ele.parameters
#     else
#         return Dict(name => ele.parameters[name] for name in names)
#     end
# end


# function get_output(ele::BaseElement, input::Dict{String,Vector{Number}})
#     # 根据input输入到模块中获取fluxes
#     fluxes = get_fluxes(ele, input)
# end