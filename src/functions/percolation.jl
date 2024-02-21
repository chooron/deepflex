function Percolation(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux(
        input_names,
        [:Percolation],
        parameters,
        percolation_func
    )
end

function percolation_func(
    input::(@NamedTuple{SoilWater::Union{T,Vector{T}}}),
    parameters::(@NamedTuple{x1::Union{T,Vector{T}}})
)::(@NamedTuple{Percolation::Union{T,Vector{T}}}) where {T<:Number}
    (Percolation=@.((parameters[:x1]^(-4)) / 4 * ((4 / 9)^(-4)) * (input[:SoilWater]^5)),)
end