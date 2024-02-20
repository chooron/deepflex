function Saturation(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux{T}(
        input_names,
        [:Saturation],
        parameters,
        saturation_func
    )
end

function saturation_func(
    input::ComponentVector{T,Vector{T},Tuple{Axis{(SoilWater=1, Rainfall=2)}}},
    parameters::ComponentVector{T,Vector{T},Tuple{Axis{(x1=1,)}}}
) where {T<:Number}
    ComponentVector(Saturation=input[:Rainfall] * (1 - (input[:SoilWater] / parameters[:x1])^2))
end