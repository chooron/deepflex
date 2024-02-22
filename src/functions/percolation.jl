function Percolation(input_names::Vector{Symbol}; parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        [:Percolation],
        parameters,
        percolation_func
    )
end

function percolation_func(
    input::gen_namedtuple_type([:SoilWater], T),
    parameters::gen_namedtuple_type([:x1], T)
)::gen_namedtuple_type([:Percolation], T) where {T<:Number}
    (Percolation=@.((parameters[:x1]^(-4)) / 4 * ((4 / 9)^(-4)) * (input[:SoilWater]^5)),)
end