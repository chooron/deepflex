struct UHFunction{uhtype}
    function UHFunction(uhtype::Symbol)
        return new{uhtype}()
    end
end

function (uh::UHFunction{:UH_1_HALF})(t, lag)
    if t - lag > 0
        typeof(lag)(1)
    else
        (t / lag)^2.5
    end
end

get_uh_tmax(::UHFunction{:UH_1_HALF}, lag) = ceil(lag)

function (uh::UHFunction{:UH_2_FULL})(t, lag)
    if t - lag * 2 > 0
        typeof(lag)(1)
    elseif t - lag > 0
        (1 - 0.5 * abs(2 - t / lag)^2.5)
    else
        (0.5 * abs(t / lag)^2.5)
    end
end

get_uh_tmax(::UHFunction{:UH_2_FULL}, lag) = 2 * ceil(lag)
