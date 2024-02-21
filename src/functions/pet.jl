function Pet(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux(
        input_names,
        [:Percolation],
        parameters,
        pet_func
    )
end

function pet_func(
    input::(@NamedTuple{Temp::Union{T,Vector{T}}, Lday::Union{T,Vector{T}}}),
    parameters::Nothing=nothing
)::(@NamedTuple{Temp::Union{T,Vector{T}}}) where {T<:Number}
    temp, lday = input[:Temp], input[:Lday]
    (Pet=@.(29.8 * lday * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)),)
end