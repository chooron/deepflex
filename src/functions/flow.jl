function Flow(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux{T}(
        input_names,
        [:Flow],
        parameters,
        flow_func
    )
end

function flow_func(
    input::ComponentVector{T,Vector{T},Tuple{Axis{(Baseflow=1, Surfaceflow=2)}}},
    parameters::Nothing=nothing
) where {T<:Number}
    ComponentVector(Flow=input[:Baseflow] + input[:Surfaceflow])
end

function flow_func(
    input::ComponentVector{T,Vector{T},Tuple{Axis{(Rainfall=1, Saturation=2, Percolation=3)}}},
    parameters::Nothing=nothing
) where {T<:Number}
    ComponentVector(Flow=input[:Rainfall] + input[:Percolation] - input[:Saturation])
end
