function baseflow(S, Smax, Qmax, f)
    @.(step_func(S) * step_func(S - Smax) * Qmax + step_func(S) * step_func(Smax - S) * Qmax * exp(-f * (Smax - S)))
end