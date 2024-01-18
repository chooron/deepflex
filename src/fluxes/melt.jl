function melt(S, Temp, Tmax, Df)
    @. step_func(Temp - Tmax) * step_func(S) * min(S, Df * (Temp - Tmax))
end