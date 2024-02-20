function Evap(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux{T}(
        input_names,
        [:Evap],
        parameters,
        evap_func
    )
end

function evap_func(
    input::ComponentVector{T,Vector{T},Tuple{Axis{(SoilWater=1, Pet=2)}}},
    parameters::ComponentVector{T,Vector{T},Tuple{Axis{(Smax=1,)}}}
) where {T<:Number}
    soil_water, pet = input[:SoilWater], input[:Pet]
    Smax = parameters[:Smax]
    ComponentVector(Evap=step_func(soil_water) * step_func(soil_water - Smax) * Pet +
                             step_func(soil_water) * step_func(Smax - soil_water) * pet * (soil_water / Smax))
end

function evap_func(
    input::ComponentVector{T,Vector{T},Tuple{Axis{(SoilWater=1, Prcp=2, Pet=3)}}},
    parameters::ComponentVector{T,Vector{T},Tuple{Axis{(x1=1,)}}}
) where {T<:Number}
    soil_water, prcp, pet = input[:SoilWater], input[:Prcp], input[:Pet]
    x1 = parameters[:x1]
    step_func(pet - prcp) * (pet - prcp) * (2 * soil_water / x1 - (soil_water / x1)^2)
end

