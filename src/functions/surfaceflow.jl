function Surfaceflow(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux(
        input_names,
        [:Surfaceflow],
        parameters,
        surfaceflow_func
    )
end

function surfaceflow_func(
    input::(@NamedTuple{SoilWater::Union{T,Vector{T}}}),
    parameters::(@NamedTuple{Smax::Union{T,Vector{T}}})
)::(@NamedTuple{Surfaceflow::Union{T,Vector{T}}}) where {T<:Number}
    soil_water = input[:SoilWater]
    Smax = parameters[:Smax]
    (Surfaceflow=(step_func(soil_water) * step_func(soil_water - Smax) * (soil_water - Smax)),)
end
