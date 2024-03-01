# Element Methods
function get_name(ele::AbstractElement)::String
    return ele.name
end

function get_output(ele::AbstractElement, input::ComponentVector)
    return nothing
end

"""
SimpleElement

SimpleElement是由多个AbstractFlux组成，无中间状态，对应RRMG的一些模型
"""
@kwdef mutable struct SimpleElement{T} <: AbstractElement where {T<:Number}
    name::String

    # parameters
    parameters::ComponentVector{T}

    # functions
    funcs::Vector{AbstractFlux}

    # attribute
    input_names::Set{Symbol}
    output_names::Set{Symbol}
end

function SimpleElement(
    ; name::String,
    parameters::ComponentVector{T},
    funcs::Vector{F},
) where {F<:AbstractFlux,T<:Number}

    input_names = Set{Symbol}()
    output_names = Set{Symbol}()

    for func in funcs
        union!(input_names, func.input_names)
        union!(output_names, func.output_names)
    end

    return SimpleElement{T}(
        name=name,
        parameters=parameters,
        funcs=funcs,
        input_names=input_names,
        output_names=output_names
    )
end

function get_output(ele::SimpleElement; input::ComponentVector{T}) where {T<:Number}
    fluxes = input
    for func in ele.funcs
        fluxes = ComponentVector(fluxes; func(fluxes)...)
    end
    return ComponentVector(; Dict(nm => fluxes[nm] for nm in ele.output_names)...)
end

"""
ODEElement

ODEElement是由多个AbstractFlux组成，有中间状态，求解方式包含continuous和discrete两种
"""
@kwdef mutable struct ODEElement{T} <: AbstractElement where {T<:Number}
    name::String

    # parameters
    parameters::ComponentVector{T}

    # states
    states::ComponentVector = ComponentVector()
    init_states::ComponentVector{T}
    state_names::Vector{Symbol}

    # functions
    funcs::Vector{AbstractFlux}

    # solve config
    get_du::Function
    solve_type::String

    # attribute
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
end

function ODEElement(
    ; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T},
    funcs::Vector{F},
    get_du::Function,
    solve_type::String=ode_solve,
) where {F<:AbstractFlux,T<:Number}

    input_names = Set{Symbol}()
    output_names = Set{Symbol}()

    for func in funcs
        if func.input_names isa Dict
            union!(input_names, Set(keys(func.input_names)))
        else
            union!(input_names, func.input_names)
        end
        union!(output_names, func.output_names)
    end

    return ODEElement{T}(
        name=name,
        parameters=parameters,
        init_states=init_states,
        state_names=collect(keys(init_states)),
        funcs=funcs,
        get_du=get_du,
        solve_type=solve_type,
        input_names=collect(input_names),
        output_names=collect(output_names)
    )
end


function get_parameters(ele::ODEElement; names::Vector{Symbol}=nothing)
    if isnothing(names)
        return ele.parameters
    else
        return Dict(name => ele.parameters[name] for name in names)
    end
end

function set_parameters!(ele::ODEElement; paraminfos::Vector{P}) where {P<:AbstractParamInfo}
    for p in paraminfos
        ele.parameters[p.name] = p.value
        for func in ele.funcs
            set_parameters!(func, paraminfos=paraminfos)
        end
    end
end

function get_states(ele::ODEElement; state_names::Union{Set{Symbol},Nothing}=nothing)
    if isnothing(state_names)
        return ele.states
    else
        available_state_names = [nm for nm in state_names if nm in ele.state_names]
        return ele.states[available_state_names]
    end
end

function pretrain!(ele::ODEElement; input::ComponentVector{T}, train_config...) where {T<:Number}
    for func in ele.funcs
        if isa(func, LuxNNFunc)
            pretrain!(func, input=input, train_config...)
        end
    end
end

function single_ele_ode_func!(du, u, p, t)
    # interpolate value by fitted functions
    tmp_ele = p[:element]
    tmp_input = ComponentVector(namedtuple(p[:input_names], [p[:itp][nm](t) for nm in p[:input_names]]))
    tmp_fluxes = get_fluxes(tmp_ele, state=u, input=tmp_input)
    tmp_du = tmp_ele.get_du(tmp_fluxes, tmp_ele.parameters)
    # return du
    for k in tmp_ele.state_names
        du[k] = tmp_du[k]
    end
end

function solve_prob!(ele::ODEElement; input::ComponentVector{T}) where {T<:Number}
    dt = 1
    cur_input_names = collect(keys(input))
    xs = 1:dt:length(input[cur_input_names[1]])
    tspan = (xs[1], xs[end])

    # fit interpolation functions
    itp_dict = Dict(nm => linear_interpolation(xs, input[nm]) for nm in cur_input_names)

    ode_parameters = (itp=itp_dict, element=ele, input_names=cur_input_names)
    # solve the problem
    prob = ODEProblem(single_ele_ode_func!, ele.init_states, tspan, ode_parameters)
    sol = solve(prob, BS3(), dt=1.0, saveat=xs, reltol=1e-3, abstol=1e-3, sensealg=ForwardDiffSensitivity())
    solved_u_matrix = hcat(sol.u...)
    solved_u = ComponentVector(namedtuple(ele.state_names, [solved_u_matrix[idx, :] for idx in 1:length(ele.state_names)]))
    ele.states = solved_u
end

function get_fluxes(ele::ODEElement; state::ComponentVector, input::ComponentVector{T}) where {T<:Number}
    fluxes = ComponentVector(input; state...)
    for func in ele.funcs
        fluxes = ComponentVector(fluxes; func(fluxes)...)
    end
    return fluxes
end

function get_output(ele::ODEElement; input::ComponentVector{T}) where {T<:Number}
    fluxes = get_fluxes(ele, state=ele.states, input=input)
    return ComponentVector(; Dict(nm => fluxes[nm] for nm in ele.output_names)...)
end


"""
LAGElement
"""
@kwdef mutable struct LAGElement{T} <: AbstractElement where {T<:Number}
    name::String

    # parameters
    parameters::ComponentVector{T}
    lag_time::ComponentVector{T}

    # states
    lag_states::ComponentVector{T}
    lag_weights::ComponentVector{T}

    # attribute
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
end

function LAGElement(
    ; name::String,
    lag_time::Union{T,ComponentVector{T}},
    lag_func::Dict{Symbol,F}
) where {T<:Number,F<:Function}

    input_names = collect(keys(lag_func))

    if typeof(lag_time) == T
        lag_time = ComponentVector(namedtuple(input_names, repeat([lag_time], length(input_names))))
    end

    parameters = ComponentVector(namedtuple(["$(k)_lag_time" for k in keys(lag_time)], [lag_time[k] for k in keys(lag_time)]))

    # init lag states
    lag_states = ComponentVector(; Dict(nm => begin
        zeros(Int(ceil(lag_time[nm])))
    end for nm in input_names)...)

    # build weight
    lag_weights = ComponentVector(; Dict(k => begin
        [
            lag_func[k](i, lag_time[k]) - lag_func[k](i - 1, lag_time[k])
            for i in 1:(ceil(lag_time[k])|>Int)
        ]
    end for k in keys(lag_time))...)

    LAGElement{T}(
        name=name,
        parameters=parameters,
        lag_time=lag_time,
        lag_states=lag_states,
        lag_weights=lag_weights,
        input_names=input_names,
        output_names=input_names)
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
    ele.lag_states = ComponentVector(; Dict(k => solved_state[k][end, :] for k in ele.input_names)...)

    output = ComponentVector(; Dict(k => solved_state[k][:, 1] for k in ele.input_names)...)
    return output
end


# function solve_prob!(ele::DCTElement; input::ComponentVector{T}) where {T<:Number}
#     data_len = length(input[first(keys(input))])
#     tmp_state = ele.init_states
#     for idx in 1:data_len
#         tmp_input = ComponentVector(; Dict(k => input[k][idx] for k in keys(input))...)
#         tmp_fluxes = get_fluxes(ele, state=tmp_state, input=tmp_input)
#         tmp_du = ele.get_du(tmp_fluxes, ele.parameters)
#         for k in keys(tmp_state)
#             update_state = tmp_state[k] + tmp_du[k]
#             tmp_state[k] = update_state
#             push!(ele.states[k], update_state)
#         end
#         tmp_state = ComponentVector(; Dict(k => tmp_state[k] + tmp_du[k] for k in keys(tmp_state))...)
#     end
# end
