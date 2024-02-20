function Pet(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux{T}(
        input_names,
        [:Percolation],
        parameters,
        pet_func
    )
end

function pet_func(
    input::ComponentVector{T,Vector{T},Tuple{Axis{(Temp=1, Lday=2)}}},
    parameters::Nothing
) where {T<:Number}
    temp, lday = input[:Temp], input[:Lday]
    ComponentVector(Pet=29.8 * lday * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2))
end