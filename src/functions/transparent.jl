function Tranparent(input_names::Vector{Symbol})
    SimpleFlux(
        input_names,
        input_names,
        nothing,
        (i::ComponentVector, p::Nothing) -> i
    )
end