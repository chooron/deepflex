function Baseflow(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux(
        input_names,
        [:Baseflow],
        parameters,
        baseflow_func
    )
end

function baseflow_func(
    input::(@NamedTuple{SoilWater::Union{T,Vector{T}}}),
    parameters::(@NamedTuple{Smax::Union{T,Vector{T}}, Qmax::Union{T,Vector{T}}, f::Union{T,Vector{T}}})
)::(@NamedTuple{Baseflow::Union{T,Vector{T}}}) where {T<:Number}
    soil_water = input[:SoilWater]
    Smax, Qmax, f = parameters[:Smax], parameters[:Qmax], parameters[:f]
    (Baseflow=@.(step_func(soil_water) * step_func(soil_water - Smax) * Qmax +
                 step_func(soil_water) * step_func(Smax - soil_water) * Qmax * exp(-f * (Smax - soil_water))),)
end

