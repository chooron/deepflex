function Smooth(conditions::Vector{T}, values::Vector{T};
    smoooth_func=step_func,
    func_args=Dict(:p1 => 5.0, :p2 => 1.0, :p3 => 0.5)) where {T<:Number}

    conditions = vcat([-Inf], conditions, [Inf])

    func = (x) -> sum([smoooth_func(x - conditions[i]; func_args...) *
                       smoooth_func(x - conditions[i+1]; func_args...) * values[i]
                       for i in 1:length(conditions)-1])
    return func
end

function step_func(x::T; p1=5.0, p2=1.0, p3=0.5) where {T<:Number}
    (tanh(p1 * x) + p2) * p3
end

