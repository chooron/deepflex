function evap(S::T, pet::T, Smax::T) where {T<:Number}
    step_func(S) * step_func(S - Smax) * pet + step_func(S) * step_func(Smax - S) * pet * (S / Smax)
end