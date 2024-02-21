function Rainfall(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux(
        input_names,
        [:Rainfall],
        parameters,
        rainfall_func
    )
end

function rainfall_func(
    input::(@NamedTuple{Prcp::Union{T,Vector{T}},Temp::Union{T,Vector{T}}}),
    parameters::(@NamedTuple{x1::Union{T,Vector{T}}})
)::(@NamedTuple{Rainfall::Union{T,Vector{T}}}) where {T<:Number}
    (Rainfall=@.(step_func(input[:Temp] - parameters[:Tmin]) * input[:Prcp]),)
end

function rainfall_func(
    input::(@NamedTuple{Prcp::Union{T,Vector{T}},Pet::Union{T,Vector{T}}}),
    parameters::Nothing
)::(@NamedTuple{Rainfall::Union{T,Vector{T}}}) where {T<:Number}
    (Rainfall=@.(step_func(input[:Prcp] - input[:Pet]) * (input[:Prcp] - input[:Pet])),)
end
