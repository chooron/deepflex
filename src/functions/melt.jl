function Melt(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    build_flux(
        input_names,
        [:Melt],
        parameters,
        melt_func
    )
end

function melt_func(
    input::gen_namedtuple_type([:SnowWater,:Temp], T),
    parameters::gen_namedtuple_type([:Tmax,:Df], T)
)::gen_namedtuple_type([:Melt], T) where {T<:Number}
    snow_water, temp = input[:SnowWater], input[:Temp]
    Tmax, Df = parameters[:Tmax], parameters[:Df]
    (Melt=@.(step_func(temp - Tmax) * step_func(snow_water) * min(snow_water, Df * (temp - Tmax))),)
end
