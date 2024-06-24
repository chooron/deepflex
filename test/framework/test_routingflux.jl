using OrdinaryDiffEq

function ifelse_func(x::Union{Num,T}) where {T<:Number}
    if x > 0.0
        return 1.0
    else
        return 0.0
    end
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

delta_t = 1.0
ts = 0:4
input = [1,2,3,4,2,1,2,2]
lag_time = 3.8

# sh, uh = uh_2_full_m(3.8, 1)
# #* 将weight作为param输入到prob中
# lag_weights = [uh_2_full(t, lag_time) for t in ts]
# lag_weights = vcat((circshift(lag_weights, -1).-lag_weights)[1:end-1])

lag_weights = [uh_1_half(t, lag_time) for t in ts]
lag_weights = vcat((circshift(lag_weights, -1).-lag_weights)[1:end-1])

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
