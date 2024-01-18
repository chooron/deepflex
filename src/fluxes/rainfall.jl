function rainfall(Prcp::T, Temp::T, Tmin::T) where{T<:Number}
    step_func(Temp - Tmin) * Prcp
end
