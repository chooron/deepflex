function ifelse_func(x::Union{Num,T}) where {T<:Number}
    if x > 0.0
        return 1.0
    else
        return 0.0
    end
end

@register_symbolic ifelse_func(x)

function step_func(x::Union{Num,T}) where {T<:Number}
    (tanh(500.0 * x) + 1.0) * 0.5
end

@register_symbolic step_func(x)