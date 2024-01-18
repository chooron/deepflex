function snowfall(Pcrp::T, Temp::T, Tmin::T) where {T<:Number}
    step_func(Tmin - Temp) * Pcrp
end