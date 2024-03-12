# Element Methods
function get_output(ele::AbstractElement, input::ComponentVector)
    return nothing
end

"""
SimpleElement

SimpleElement是由多个AbstractFlux组成，无中间状态，对应RRMG的一些模型
"""
struct SimpleElement <: AbstractElement
    name::Symbol

    # attribute
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
    parameters_names::Vector{Symbol}

    # functions
    funcs::Vector{AbstractFlux}
end

function SimpleElement(
    ; name::Symbol,
    funcs::Vector{F},
) where {F<:AbstractFlux}

    input_names, output_names, parameter_names = get_func_infos(funcs)

    parameters = ComponentVector{Number}(namedtuple(parameter_names, ones(Number, length(parameter_names))))

    return SimpleElement(
        name,
        input_names,
        output_names,
        parameters_names,
        funcs,
        parameters
    )
end

function get_output(ele::SimpleElement; input::ComponentVector{T}, parameters::ComponentVector{T}) where {T<:Number}
    fluxes = input
    for func in ele.funcs
        fluxes = ComponentVector(fluxes; func(fluxes, parameters)...)
    end
    return fluxes[ele.output_names]
end

"""
ODEElement

ODEElement是由多个AbstractFlux组成，有中间状态，求解方式包含continuous和discrete两种
"""
struct ODEElement <: AbstractElement
    name::Symbol

    # attribute
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
    parameter_names::Vector{Symbol}
    state_names::Vector{Symbol}

    # functions
    funcs::Vector{AbstractFlux}
    d_funcs::Vector{AbstractFlux}
end

function ODEElement(
    name::Symbol;
    funcs::Vector{F},
    d_funcs::Vector{F},
) where {F<:AbstractFlux}

    input_names1, output_names, parameter_names1 = get_func_infos(funcs)
    input_names2, states_names, parameter_names2 = get_d_func_infos(d_funcs)

    input_names = union(input_names1, input_names2)
    parameter_names = union(parameter_names1, parameter_names2)

    return ODEElement(
        name,
        input_names,
        output_names,
        parameter_names,
        states_names,
        funcs,
        d_funcs
    )
end

function set_solver!(ele::ODEElement, solver::AbstractSolver)
    setproperty!(ele, :solver, solver)
end

function pretrain!(ele::ODEElement; input::ComponentVector{T}, train_config...) where {T<:Number}
    for func in ele.funcs
        if isa(func, LuxNNFunc)
            pretrain!(func, input=input, train_config...)
        end
    end
end

function solve_prob(
    ele::ODEElement;
    input::ComponentVector{T},
    parameters::ComponentVector{T},
    init_states::ComponentVector{T},
    solver::AbstractSolver,
)::ComponentVector{T} where {T<:Number}
    # fit interpolation functions
    itp_dict = Dict(nm => linear_interpolation(input[:time], input[nm]) for nm in element.input_names)
    # solve the problem
    function singel_ele_ode_func!(du, u, p)
        tmp_input = ComponentVector(namedtuple(ele.input_names, [itp_dict[nm](t) for nm in element.input_names]))
        tmp_fluxes = get_output(ele, input=tmp_input, state=u, parameters=p)
        for d_func in tmp_ele.d_funcs
            du[d_func.output_names[1]] = d_func(tmp_fluxes)[d_func.output_names[1]]
        end
    end
    # return solved result
    solver(singel_ele_ode_func!, parameters, init_states, time_config)
end

function get_output(
    ele::ODEElement;
    input::ComponentVector{T},
    states::ComponentVector{T},
    parameters::ComponentVector{T}
)::ComponentVector{T} where {T<:Number}
    fluxes = ComponentVector(input; states...)
    for func in ele.funcs
        fluxes = ComponentVector(fluxes; func(fluxes, parameters)...)
    end
    return fluxes
end

"""
LAGElement
"""
@kwdef mutable struct LAGElement <: AbstractElement
    name::Symbol

    # attribute
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}

    # func
    lag_func::Dict{Symbol,F}
    step_func::Function

    # parameters
    lag_time::ComponentVector

    # lag states
    lag_states::ComponentVector

    # lag weights
    lag_weights::ComponentVector
end

function LAGElement(
    name::Symbol;
    lag_func::Dict{Symbol,F},
    step_func::Function=DEFAULT_SMOOTHER
) where {F<:Function}

    input_names = Set(keys(lag_func))

    LAGElement(
        name,
        input_names,
        input_names,
        lag_func,
        step_func,
        ComponentVector(),
        ComponentVector(),
        ComponentVector()
    )
end

function set_lag_time!(ele::LAGElement; lag_time::ComponentVector{T}) where {T<:Number}
    if Set(keys(lag_time)) != ele.input_names
        @error "$(Set(key(lag_time))) is not consistent with the states of element($(ele.name)): $(ele.input_names)"
    end
    setfield!(ele, :lag_time, lag_time)
    # init lag states
    lag_states = ComponentVector(namedtuple(input_names,
        [zeros(Int(ceil(lag_time[nm]))) for nm in ele.input_names]))

    # build weight
    lag_weights = ComponentVector(namedtuple(input_names,
        [[lag_func[k](T(i), T(lag_time[k]), step_func) -
          lag_func[k](T(i - 1), T(lag_time[k]), step_func)
          for i in 1:(ceil(lag_time[k])|>Int)]
         for k in ele.input_names]))

    setfield!(ele, :lag_states, lag_states)
    setfield!(ele, :lag_weights, lag_weights)
end

function solve_lag(ele::LAGElement; input::ComponentVector{T}) where {T<:Number}
    max_weight_len = maximum([length(ele.lag_weights[k]) for k in keys(ele.lag_weights)])
    max_input_len = maximum([length(input[k]) for k in keys(input)])

    output = Dict(k => zeros(T, max_input_len, max_weight_len) for k in keys(input))

    for k in keys(output)
        w, ls, ip = ele.lag_weights[k], ele.lag_states[k], input[k]
        for t in 1:max_input_len
            updated_state = ls .+ ip[t] .* w
            output[k][t, 1:length(w)] .= updated_state
            ls = vcat(updated_state[2:end], 0)
        end
    end
    return output
end

function get_output(ele::LAGElement; input::ComponentVector{T}) where {T<:Number}
    input = input[collect(ele.input_names)]
    solved_state = solve_lag(ele, input=input)

    for k in ele.input_names
        solved_state[k][end, 1:end-1] = solved_state[k][end, 2:end]
        solved_state[k][end, end] = 0
        solved_state[k][end, :]
    end

    # Get the new lag value to restart
    setfield!(ele, :lag_states, ComponentVector(namedtuple(ele.input_names, [solved_state[k][end, :] for k in ele.input_names])))
    ComponentVector(namedtuple(ele.input_names, [solved_state[k][:, 1] for k in ele.input_names]))
end