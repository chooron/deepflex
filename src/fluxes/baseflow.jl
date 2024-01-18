function baseflow(S::T, Smax::T, Qmax::T, f::T) where {T<:Number}
    step_func(S) * step_func(S - Smax) * Qmax + step_func(S) * step_func(Smax - S) * Qmax * exp(-f * (Smax - S))
end