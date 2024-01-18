function surfaceflow(S, Smax)
    step_func(S) * step_func(S .- Smax) .* (S .- Smax)
end