function Rainfall(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux{T}(
        input_names,
        [:Rainfall],
        parameters,
        rainfall_func
    )
end

function rainfall_func(
    input::ComponentVector{T,Vector{T},Tuple{Axis{(Prcp=1, Temp=2)}}},
    parameters::ComponentVector{T,Vector{T},Tuple{Axis{(Tmin=1,)}}}
) where {T<:Number}
    ComponentVector(Rainfall=step_func(input[:Temp] - parameters[:Tmin]) * input[:Prcp])
end

function rainfall_func(
    input::ComponentVector{T,Vector{T},Tuple{Axis{(Prcp=1, Pet=2)}}},
    parameters::Nothing
) where {T<:Number}
    ComponentVector(Rainfall=step_func(input[:Prcp] - input[:Pet]) * (input[:Prcp] - input[:Pet]))
end
