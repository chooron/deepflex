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