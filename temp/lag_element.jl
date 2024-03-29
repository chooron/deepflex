
"""
LAGElement
"""
struct LAGElement <: AbstractElement
    name::Symbol

    # attribute
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}

    # parameters names
    param_names::Vector{Symbol}
    lagtime_dict::Dict{Symbol,Symbol}

    # func
    lag_func::Dict{Symbol,Function}
    step_func::Function
end

function LAGElement(
    ; name::Symbol,
    lagtime_dict::Dict{Symbol,Symbol},
    lag_func::Dict{Symbol,F}
) where {F<:Function}

    step_func = DEFAULT_SMOOTHER

    param_names = Vector{Symbol}()
    for v in values(lagtime_dict)
        union!(param_names, [v])
    end

    input_names = collect(keys(lag_func))

    LAGElement(
        name,
        input_names,
        input_names,
        param_names,
        lagtime_dict,
        lag_func,
        step_func,
    )
end

function preprocess_parameters(ele::LAGElement; parameters::NamedTuple)
    # init lag states
    lag_states = namedtuple(ele.input_names, [zeros(Int(ceil(parameters[ele.lagtime_dict[nm]]))) for nm in ele.input_names])

    # build weight
    lag_weights = namedtuple(ele.input_names,
        [[ele.lag_func[nm](T(i), T(parameters[ele.lagtime_dict[nm]]), ele.step_func) -
          ele.lag_func[nm](T(i - 1), T(parameters[ele.lagtime_dict[nm]]), ele.step_func)
          for i in 1:(ceil(parameters[ele.lagtime_dict[nm]])|>Int)] for nm in ele.input_names])

    lag_states, lag_weights
end

function solve_lag(input::NamedTuple, lag_states::NamedTuple, lag_weights::NamedTuple)
    max_weight_len = maximum([length(lag_weights[k]) for k in keys(lag_weights)])
    max_input_len = maximum([length(input[k]) for k in keys(input)])

    output = Dict(k => zeros(max_input_len, max_weight_len) for k in keys(input))

    for k in keys(output)
        w, ls, ip = lag_weights[k], lag_states[k], input[k]
        for t in 1:max_input_len
            updated_state = ls .+ ip[t] .* w
            output[k][t, 1:length(w)] .= updated_state
            ls = vcat(updated_state[2:end], 0)
        end
    end
    return output
end

function get_output(ele::LAGElement; input::NamedTuple, parameters::NamedTuple)
    lag_states, lag_weights = preprocess_parameters(ele, parameters=parameters)

    input = input[ele.input_names]
    solved_state = solve_lag(input, lag_states, lag_weights)

    for k in ele.input_names
        solved_state[k][end, 1:end-1] = solved_state[k][end, 2:end]
        solved_state[k][end, end] = 0
        solved_state[k][end, :]
    end
    # lag element don't need save states
    # updated_lag_states = NamedTuple(namedtuple(ele.input_names, [solved_state[k][end, :] for k in ele.input_names]))
    namedtuple(ele.input_names, [solved_state[k][:, 1] for k in ele.input_names])
end