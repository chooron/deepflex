function snowfall(; P::Vector{float}, T::Vector{float}, Tmin::float)
    step_fct(Tmin - T) * P
end