function Baseflow(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux{T}(
        input_names,
        [:Baseflow],
        parameters,
        baseflow_func
    )
end

function baseflow_func(
    input::ComponentVector{T,Vector{T},Tuple{Axis{(SoilWater=1,)}}},
    parameters::ComponentVector{T,Vector{T},Tuple{Axis{(Smax=1, Qmax=2, f=3)}}}
) where {T<:Number}
    soil_water = input[:SoilWater]
    Smax, Qmax, f = parameters[:Smax], parameters[:Qmax], parameters[:f]
    ComponentVector(Baseflow=step_func(soil_water) * step_func(soil_water - Smax) * Qmax +
                             step_func(soil_water) * step_func(Smax - soil_water) * Qmax * exp(-f * (Smax - soil_water)))
end

