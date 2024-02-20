function Melt(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux{T}(
        input_names,
        [:Melt],
        parameters,
        melt_func
    )
end

function melt_func(
    input::ComponentVector{T,Vector{T},Tuple{Axis{(SnowWater=1, Temp=2)}}},
    parameters::ComponentVector{T,Vector{T},Tuple{Axis{(Tmax=1, Df=2)}}}
) where {T<:Number}
    snow_water, temp = input[:SnowWater], input[:Temp]
    Tmax, Df = parameters[:Tmax], parameters[:Df]
    ComponentVector(Melt=step_func(temp - Tmax) * step_func(snow_water) * min(snow_water, Df * (temp - Tmax)))
end
