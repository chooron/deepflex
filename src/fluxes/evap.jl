function evap(S, pet, Smax)
    (step_fct(S) * step_fct(S - Smax) * pet + step_fct(S) * step_fct(Smax - S) * pet * (S / Smax))
end