function Saturation(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux(
        input_names,
        [:Saturation],
        parameters,
        saturation_func
    )
end

function saturation_func(
    input::(@NamedTuple{SoilWater::Union{T,Vector{T}},Rainfall::Union{T,Vector{T}}}),
    parameters::(@NamedTuple{x1::Union{T,Vector{T}}})
)::(@NamedTuple{Saturation::Union{T,Vector{T}}}) where {T<:Number}
    (Saturation=@.(input[:Rainfall] * (1 - (input[:SoilWater] / parameters[:x1])^2)),)
end