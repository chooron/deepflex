function Rainfall(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    build_flux(
        input_names,
        [:Rainfall],
        parameters,
        rainfall_func
    )
end

function rainfall_func(
    input::gen_namedtuple_type([:Prcp,:Temp], T),
    parameters::gen_namedtuple_type([:Tmin], T)
)::gen_namedtuple_type([:Rainfall], T) where {T<:Number}
    (Rainfall=@.(step_func(input[:Temp] - parameters[:Tmin]) * input[:Prcp]),)
end

function rainfall_func(
    input::gen_namedtuple_type([:Prcp,:Pet], T),
    parameters::Nothing=nothing
)::gen_namedtuple_type([:Rainfall], T) where {T<:Number}
    (Rainfall=@.(step_func(input[:Prcp] - input[:Pet]) * (input[:Prcp] - input[:Pet])),)
end
