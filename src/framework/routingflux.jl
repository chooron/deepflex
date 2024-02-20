mutable struct RoutingFlux{T<:Number} <: AbstractFlux

    # attribute
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}

    # parameters
    lag_time::ComponentVector{T}

    # states
    lag_states::ComponentVector{T}

    # functions
    lag_func::Function

    # weight
    lag_weights::ComponentVector{T}

end

function RoutingFlux(input_names::Vector{Symbol}; lag_time::Union{T,ComponentVector{T}}, lag_func::Function) where {T<:Number}
    if typeof(lag_time) == T
        lag_time = ComponentVector(; Dict(nm => lag_time for nm in input_names)...)
    end

    # init lag states
    lag_states = ComponentVector(; Dict(nm => begin
        zeros(Int(ceil(lag_time[nm])))
    end for nm in input_names)...)
    # build weight
    lag_weights = ComponentVector(; Dict(k => begin
        [
            lag_func(i + 1, ceil(lag_time[k])) - lag_func(i, ceil(lag_time[k]))
            for i in 1:(ceil(lag_time[k])|>Int)
        ]
    end for k in keys(lag_time))...)

    RoutingFlux{T}(
        input_names,
        input_names,
        lag_time,
        lag_states,
        lag_func,
        lag_weights)
end

function solve_lag(flux::RoutingFlux; input::ComponentVector{T}, lag_state::ComponentVector{T}) where {T<:Number}
    max_weight_len = max([length(weight(k)) for k in keys(flux.weight)])
    max_input_len = max([length(input(k)) for k in keys(input)])
    output = ComponentVector(; Dict(k => zeros(T, max_input_len, max_weight_len))...)

    for k in keys(output)
        for (w, ls, i) in zip(flux.weight[k], lag_state[k], input[k])
            for ts in 1:max_input_len
                updated_state = ls .+ i[ts] .* w
                output[k][ts, 1:length(w)] = updated_state
                ls = vcat(updated_state[2:end], 0)
            end
        end
    end
    return output
end

function get_output(flux::RoutingFlux; input::ComponentVector{T}) where {T<:Number}
    solved_state = solve_lag(flux, input=input, lag_state=flux.lag_state)

    # Get the new lag value to restart
    final_states = ComponentVector(; Dict(k => begin
        tmp_state = solved_state[k][end, :]
        tmp_state[:, 1:end-1] = tmp_state[:, 2:end]
        tmp_state[:, end] = T(0)
    end for k in keys(flux.input_names))...)

    flux.lag_states = final_states

    output = ComponentVector(; Dict(k => solved_state[k][:, 1] for k in flux.input_names)...)
    return output
end


function unit_hydro1(bin, len)
    value = begin
        step_func(bin - len) +
        step_func(len - bin) * step_func(bin) * (bin / len)^2.5
    end
    return value
end

function unit_hydro2(bin, len)
    half_len = len / 2
    value = begin
        step_func(bin - len) * 1 +
        step_func(len - bin) * step_func(bin - half_len) * (1 - 0.5 * abs(2 - bin / half_len)^2.5) +
        step_func(half_len - bin) * step_func(bin) * (0.5 * abs(bin / half_len)^2.5)
    end
    return value
end