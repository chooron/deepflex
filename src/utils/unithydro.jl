function init_lag_state(lag_func::Function, step_func::Function, lag_time::T, delta_t::T) where {T<:Number}
    delay = lag_time / delta_t
    ts = 1:(ceil(delay)|>Int)
    lag_state = zeros(Num, (2, length(ts)))
    tmp_state = [lag_func(t, lag_time, step_func) for t in ts]
    lag_state[1, :] = vcat([tmp_state[1]], (circshift(tmp_state, -1).-tmp_state)[1:end-1])
    lag_state
end

function update_lag_state!(lag_state::Array{T}, input::Number) where {T<:Number}
    lag_state[2, :] = lag_state[1, :] .* input + lag_state[2, :]
    lag_state[2, :] = circshift(lag_state[2, :], -1)
    lag_state[2, end] = 0
end

function solve_lag_flux(input::Vector, lag_time::Number, lag_func::Function; kwargs...)
    delta_t = 1.0
    # ts = 1:(ceil(lag_time / delta_t)|>Int)
    ts = 0:200
    #* 将weight作为param输入到prob中
    lag_weights = [lag_func(t, lag_time) for t in ts]
    lag_weights = vcat([lag_weights[1]], (circshift(lag_weights, -1).-lag_weights)[1:end-1])

    #* 首先将lagflux转换为discrete problem
    function discrete_prob(u, p, t)
        u = circshift(u, -1)
        u[end] = 0.0
        tmp_u = input[Int(t)] .* p .+ u
        tmp_u
    end

    prob = DiscreteProblem(discrete_prob, lag_weights, (1.0, length(input)), lag_weights)
    #* 求解这个问题
    sol = solve(prob, FunctionMap())
    #* 得到权重计算结果
    sol[1, :]
end

function uh_1_half(tdx, lag; kw...)
    sf = ifelse_func
    value = @.(sf(tdx - lag) +
               sf(lag - tdx) * sf(tdx) * (tdx / lag)^2.5)
    return value
end

function uh_2_full(tdx, lag; kw...)
    sf = ifelse_func
    double_lag = lag * 2
    
    value = @.(sf(tdx - double_lag) * 1 +
               sf(double_lag - tdx) * sf(tdx - lag) * (1 - 0.5 * abs(2 - tdx / lag)^2.5) +
               sf(lag - tdx) * (0.5 * abs(tdx / lag)^2.5))
    return value
end

function uh_3_half(input, lag; kw...)
    sf = get(kw, :smooth_func, step_func)
    ff = @.(1 / (0.5 * delay^2))
    value = @.(sf(lag - input) * ff * (0.5 * input^2 - 0.5 * (input - 1)^2) +
               sf(input - lag) * (0.5 * delay^2 - 0.5 * (t - 1)^2))
    return value
end

function uh_4_full(input, lag; kw...)
    ff = @.(0.5 / (0.5 * (0.5 * lag)^2))
    half_lag = 0.5 * lag
    max(ff .* (input - half_lag) .* sign(half_lag - input) + ff .* half_lag, 0)
end

function uh_5_half(input, lag; kw...)
    stepsize = Int32(7 ÷ lag)
    limits = 0:stepsize:7
    limits[end] = 7
    exp(-input)
end