function uh_1_half(bin, len, sf)
    value = begin
        sf(bin - len) +
        sf(len - bin) * sf(bin) * (bin / len)^2.5
    end
    return value
end

function uh_2_full(bin, len, sf)
    half_len = len / 2
    value = begin
        sf(bin - len) * 1 +
        sf(len - bin) * sf(bin - half_len) * (1 - 0.5 * abs(2 - bin / half_len)^2.5) +
        sf(half_len - bin) * sf(bin) * (0.5 * abs(bin / half_len)^2.5)
    end
    return value
end

function uh_3_half(bin, len, sf)
    ff = @.(1 / (0.5 * delay^2))
    value = begin
        @.(sf(len - bin) * ff * (0.5 * bin^2 - 0.5 * (bin - 1)^2) +
           sf(bin - len) * (0.5 * delay^2 - 0.5 * (t - 1)^2))
    end
    return value
end

function uh_4_full(bin, len, sf)
    ff = @.(0.5 / (0.5 * (0.5 * len)^2))
    half_len = 0.5 * len
    max(ff .* (bin - half_len) .* sign(half_len - bin) + ff .* half_len, 0)
end

function uh_5_half(bin, len, sf)
    stepsize = Int32(7 ÷ len)
    limits = 0:stepsize:7
    limits[end] = 7
    exp(-bin)
end