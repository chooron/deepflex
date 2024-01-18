function rainfall(Prcp::Vector{T}, Temp::Vector{T}, Tmin::T) where{T<:Number}
    @.step_func(Temp - Tmin) * Prcp
end
