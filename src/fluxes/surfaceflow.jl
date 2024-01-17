function surfaceflow(S, Smax)
    step_fct(S) * step_fct(S - Smax) * (S - Smax)
end