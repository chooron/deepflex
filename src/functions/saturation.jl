function Saturation(input_names::Vector{Symbol}; parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        [:Saturation],
        parameters,
        saturation_func
    )
end

function saturation_func(
    input::gen_namedtuple_type([:SoilWater,:Rainfall], T),
    parameters::gen_namedtuple_type([:x1], T)
)::gen_namedtuple_type([:Saturation], T) where {T<:Number}
    (Saturation=@.(input[:Rainfall] * (1 - (input[:SoilWater] / parameters[:x1])^2)),)
end