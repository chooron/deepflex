function melt(S, T; Tmax, Df)
    step_fct(T - Tmax) * step_fct(S) * minimum(S, Df * (T - Tmax))
end