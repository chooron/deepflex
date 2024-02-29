function Percolation(input_names::Union{Vector{Symbol},Vector{Dict{Symbol,Symbol}}},
    output_names::Vector{Symbol}=[:Percolation];
    parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        output_names,
        parameters=parameters,
        func=percolation_func
    )
end

function percolation_func(
    input::gen_namedtuple_type([:SoilWater], T),
    parameters::gen_namedtuple_type([:x1], T)
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    [@.((parameters[:x1]^(-4)) / 4 * ((4 / 9)^(-4)) * (input[:SoilWater]^5))]
end