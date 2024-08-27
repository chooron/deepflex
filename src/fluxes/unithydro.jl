function uh_1_half(lag; kw...)
    timeidx = 1:ceil(lag)
    sf = get(kw, :smooth_func, ifelse_func)
    value = @.(sf(timeidx - lag) +
               sf(lag - timeidx) * sf(timeidx) * (timeidx / lag)^2.5)
    vcat([value[1]], value[2:end] .- value[1:end-1])
end

function uh_2_full(lag; kw...)
    sf = get(kw, :smooth_func, ifelse_func)
    double_lag = lag * 2
    timeidx = 1:ceil(double_lag)

    value = @.(sf(timeidx - double_lag) * 1 +
               sf(double_lag - timeidx) * sf(timeidx - lag) * (1 - 0.5 * abs(2 - timeidx / lag)^2.5) +
               sf(lag - timeidx) * (0.5 * abs(timeidx / lag)^2.5))

    vcat(value[1], value[2:end] .- value[1:end-1])
end

function uh_3_half(lag; kw...)
    sf = get(kw, :smooth_func, ifelse_func)

    timeidx = 1:ceil(lag)
    ff = 1 / (0.5 * lag^2)
    value = @.(sf(lag - timeidx) * ff * (0.5 * timeidx^2 - 0.5 * (timeidx - 1)^2) +
               sf(timeidx - lag) * ff * (0.5 * lag^2 - 0.5 * (timeidx - 1)^2))
    return value
end

