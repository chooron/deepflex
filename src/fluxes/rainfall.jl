function rainfall(; P::Vector{float}, T::Vector{float}, Tmin::Vector{float})
    step_fct(T - Tmin) * P
end
