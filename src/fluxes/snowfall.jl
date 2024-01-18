function snowfall(Pcrp::Vector{T}, Temp::Vector{T}, Tmin::T) where {T<:Number}
    @.step_func(Tmin - Temp) * Pcrp
end