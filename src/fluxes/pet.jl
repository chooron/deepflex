function pet(Temp::T, Lday::T) where {T<:Number}
    29.8 * Lday * 0.611 * exp((17.3 * Temp) / (Temp + 237.3)) / (Temp + 273.2)
end