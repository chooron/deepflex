function Snowfall(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux{T}(
        input_names,
        [:Snowfall],
        parameters,
        snowfall_func
    )
end

function snowfall_func(
    input::ComponentVector{T,Vector{T},Tuple{Axis{(Prcp=1, Temp=2)}}},
    parameters::ComponentVector{T,Vector{T},Tuple{Axis{(Tmin=1,)}}}
) where {T<:Number}
    ComponentVector(Snowfall=step_func(parameters[:Tmin] - input[:Temp]) * input[:Prcp])
end
