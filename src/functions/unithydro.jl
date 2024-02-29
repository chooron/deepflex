function uh_1_half(bin, len)
    value = begin
        step_func(bin - len) +
        step_func(len - bin) * step_func(bin) * (bin / len)^2.5
    end
    return value
end

function uh_2_full(bin, len)
    half_len = len / 2
    value = begin
        step_func(bin - len) * 1 +
        step_func(len - bin) * step_func(bin - half_len) * (1 - 0.5 * abs(2 - bin / half_len)^2.5) +
        step_func(half_len - bin) * step_func(bin) * (0.5 * abs(bin / half_len)^2.5)
    end
    return value
end

function uh_3_half(bin, len)
    ff = @.(1 / (0.5 * delay^2))
    value = begin
        @.(step_func(len - bin) * ff * (0.5 * bin^2 - 0.5 * (bin - 1)^2) +
           step_func(bin - len) * (0.5 * delay^2 - 0.5 * (t - 1)^2))
    end
    return value
end

function uh_4_full(bin, len)
    ff = @.(0.5 / (0.5 * (0.5 * len)^2))
    half_len = 0.5 * len
    max(ff .* (bin - half_len) .* sign(half_len - bin) + ff .* half_len, 0)
end

function uh_5_half(bin, len)
    stepsize = Int32(7 ÷ len)
    limits = 0:stepsize:7
    limits[end] = 7
    exp(-bin)
end