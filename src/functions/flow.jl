function Flow(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux(
        input_names,
        [:Flow],
        parameters,
        flow_func
    )
end

function flow_func(
    input::(@NamedTuple{Baseflow::Union{T,Vector{T}}, Surfaceflow::Union{T,Vector{T}}}),
    parameters::Nothing=nothing
)::(@NamedTuple{Flow::Union{T,Vector{T}}}) where {T<:Number}
    (Flow=input[:Baseflow] .+ input[:Surfaceflow],)
end

function flow_func(
    input::(@NamedTuple{Rainfall::Union{T,Vector{T}}, Percolation::Union{T,Vector{T}}, Saturation::Union{T,Vector{T}}}),
    parameters::Nothing=nothing
)::(@NamedTuple{Flow::Union{T,Vector{T}}}) where {T<:Number}
    (Flow=input[:Rainfall] .+ input[:Percolation] .- input[:Saturation],)
end
