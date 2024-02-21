function Snowfall(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux(
        input_names,
        [:Snowfall],
        parameters,
        snowfall_func
    )
end

function snowfall_func(
    input::(@NamedTuple{Prcp::Union{T,Vector{T}},Temp::Union{T,Vector{T}}}),
    parameters::(@NamedTuple{Tmin::Union{T,Vector{T}}})
)::(@NamedTuple{Snowfall::Union{T,Vector{T}}}) where {T<:Number}
    (Snowfall=@.(step_func(parameters[:Tmin] - input[:Temp]) * input[:Prcp]),)
end
