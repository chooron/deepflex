function melt(S::T, Temp::T, Tmax::T, Df::T) where {T<:Number}
    step_func(Temp - Tmax) * step_func(S) * min(S, Df * (Temp - Tmax))
end