function Evap(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux(
        input_names,
        [:Evap],
        parameters,
        evap_func
    )
end

function evap_func(
    input::(@NamedTuple{SoilWater::Union{T,Vector{T}}, Pet::Union{T,Vector{T}}}),
    parameters::(@NamedTuple{Smax::Union{T,Vector{T}}})
)::(@NamedTuple{Evap::Union{T,Vector{T}}}) where {T<:Number}
    soil_water, pet = input[:SoilWater], input[:Pet]
    Smax = parameters[:Smax]
    (Evap=@.(step_func(soil_water) * step_func(soil_water - Smax) * pet +
             step_func(soil_water) * step_func(Smax - soil_water) * pet * (soil_water / Smax)),)
end

function evap_func(
    input::(@NamedTuple{SoilWater::Union{T,Vector{T}}, Prcp::Union{T,Vector{T}}, Pet::Union{T,Vector{T}}}),
    parameters::(@NamedTuple{x1::Union{T,Vector{T}}})
)::NamedTuple{Evap::Union{T,Vector{T}}} where {T<:Number}
    soil_water, prcp, pet = input[:SoilWater], input[:Prcp], input[:Pet]
    x1 = parameters[:x1]
    (Evap=@.(step_func(pet - prcp) * (pet - prcp) * (2 * soil_water / x1 - (soil_water / x1)^2)),)
end

