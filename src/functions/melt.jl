function Melt(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux(
        input_names,
        [:Melt],
        parameters,
        melt_func
    )
end

function melt_func(
    input::(@NamedTuple{SnowWater::Union{T,Vector{T}}, Temp::Union{T,Vector{T}}}),
    parameters::(@NamedTuple{Tmax::Union{T,Vector{T}}, Df::Union{T,Vector{T}}})
)::(@NamedTuple{Melt::Union{T,Vector{T}}, Df::Union{T,Vector{T}}}) where {T<:Number}
    snow_water, temp = input[:SnowWater], input[:Temp]
    Tmax, Df = parameters[:Tmax], parameters[:Df]
    (Melt=@.(step_func(temp - Tmax) * step_func(snow_water) * min(snow_water, Df * (temp - Tmax))),)
end
