function Pet(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    build_flux(
        input_names,
        [:Pet],
        parameters,
        pet_func
    )
end

function pet_func(
    input::gen_namedtuple_type([:Temp,:Lday], T),
    parameters::Nothing=nothing,
)::gen_namedtuple_type([:Pet], T) where {T<:Number}
    temp, lday = input[:Temp], input[:Lday]
    (Pet=@.(29.8 * lday * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)),)
end