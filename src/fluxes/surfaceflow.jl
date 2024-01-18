function surfaceflow(S::T, Smax::T) where {T<:Number}
    step_func(S) * step_func(S .- Smax) .* (S .- Smax)
end