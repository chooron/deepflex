function Tranparent(input_names::Vector{Symbol}; parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        input_names,
        parameters,
        (i::NamedTuple, parameters::Nothing=nothing) -> i
    )
end