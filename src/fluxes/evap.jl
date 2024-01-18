function evap(S, pet, Smax)
    @. (step_func(S) * step_func(S - Smax) * pet + step_func(S) * step_func(Smax - S) * pet * (S / Smax))
end