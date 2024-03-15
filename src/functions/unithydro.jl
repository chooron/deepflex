function init_lag_state(lag_func::Function, lag_time::T, delta_t::T) where {T<:Number}
    delay = lag_time / delta_t
    lag_state = zeros(T, (2, length(tt)))
    lag_state[1, :] = vcat(T[0], [lag_func(t + 1, lag_time) - lag_func(t, lag_time) for t in 1:(ceil(delay)|>Int)])
    lag_state
end

function update_lag_state!(lag_state::Vector{T}, input::T) where {T<:Number}
    lag_state[2, :] = lag_state[1, :] .* input + lag_state[2, :]
    circshift!(lag_state[2, :], -1)
    lag_state[2, end] = 0
end

function uh_1_half(input, lag, sf)
    value = @.(sf(input - lag) +
               sf(lag - input) * sf(input) * (input / lag)^2.5)
    return value
end

function uh_2_full(input, lag, sf)
    half_lag = lag ./ 2
    value = @.(sf(input - lag) * 1 +
               sf(lag - input) * sf(input - half_lag) * (1 - 0.5 * abs(2 - input / half_lag)^2.5) +
               sf(half_lag - input) * sf(input) * (0.5 * abs(input / half_lag)^2.5))
    return value
end

function uh_3_half(input, lag, sf)
    ff = @.(1 / (0.5 * delay^2))
    value = @.(sf(lag - input) * ff * (0.5 * input^2 - 0.5 * (input - 1)^2) +
               sf(input - lag) * (0.5 * delay^2 - 0.5 * (t - 1)^2))
    return value
end

function uh_4_full(input, lag, sf)
    ff = @.(0.5 / (0.5 * (0.5 * lag)^2))
    half_lag = 0.5 * lag
    max(ff .* (input - half_lag) .* sign(half_lag - input) + ff .* half_lag, 0)
end

function uh_5_half(input, lag, sf)
    stepsize = Int32(7 ÷ lag)
    limits = 0:stepsize:7
    limits[end] = 7
    exp(-input)
end