function Surfaceflow(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux{T}(
        input_names,
        [:Surfaceflow],
        parameters,
        surfaceflow_func
    )
end

function surfaceflow_func(
    input::ComponentVector{T,Vector{T},Tuple{Axis{(SoilWater=1,)}}},
    parameters::ComponentVector{T,Vector{T},Tuple{Axis{(Smax=1,)}}}
) where {T<:Number}
    soil_water = input[:SoilWater]
    Smax = parameters[:Smax]
    ComponentVector(Surfaceflow=step_func(soil_water) * step_func(soil_water - Smax) * (soil_water - Smax))
end