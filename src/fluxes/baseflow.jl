function baseflow(S, Smax, Qmax, f)
    (step_fct(S) * step_fct(S - Smax) * Qmax + step_fct(S) * step_fct(Smax - S) * Qmax * exp(-f * (Smax - S)))
end