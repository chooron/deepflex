function Percolation(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux{T}(
        input_names,
        [:Percolation],
        parameters,
        percolation_func
    )
end

function percolation_func(
    input::ComponentVector{T,Vector{T},Tuple{Axis{(SoilWater=1,)}}},
    parameters::ComponentVector{T,Vector{T},Tuple{Axis{(x1=1,)}}}
) where {T<:Number}
    ComponentVector(Percolation=(parameters[:x1]^(-4)) / 4 * ((4 / 9)^(-4)) * (input[:SoilWater]^5))
end