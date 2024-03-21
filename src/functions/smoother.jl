function step_function(parameter::NamedTuple)
    function inner_step_func(x::T) where {T<:Number}
        _step_func(x, parameter)
    end
    function inner_step_func(x::Vector{T}) where {T<:Number}
        _step_func.(x, [parameter])
    end
    inner_step_func
end

function _step_func(x::Union{Num,T}, p::(@NamedTuple{v::T})) where {T<:Number}
    if x > p.v
        return 1
    else
        return 0
    end
end

function _step_func(x::Union{Num,T}, p::(@NamedTuple{p1::T, p2::T, p3::T})) where {T<:Number}
    (tanh(p.p1 * x) + p.p2) * p.p3
end

const DEFAULT_SMOOTHER = step_function((p1=5.0, p2=1.0, p3=0.5))
const DEFAULT_NO_SMOOTHER = step_function((v=0.0,))
