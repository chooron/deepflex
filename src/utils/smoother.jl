# function step_function(parameter::NamedTuple)
#     function inner_step_func(x::T) where {T<:Number}
#         _step_func(x, parameter)
#     end
#     function inner_step_func(x::Vector{T}) where {T<:Number}
#         _step_func.(x, [parameter])
#     end
#     inner_step_func
# end

function ifelse_func(x::Union{Num,T}) where {T<:Number}
    if x > 0.0
        return 1.0
    else
        return 0.0
    end
end

@register_symbolic ifelse_func(x)

function step_func(x::Union{Num,T}) where {T<:Number}
    (tanh(5.0 * x) + 1.0) * 0.5
end

@register_symbolic step_func(x)