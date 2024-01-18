function step_func(x::T; p1=5.0, p2=1.0, p3=0.5) where {T<:Number}
    (tanh(p1 * x) + p2) * p3
end

