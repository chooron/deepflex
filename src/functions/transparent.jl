function Tranparent(input_names::Vector{Symbol}; parameters::Nothing=nothing)
    build_flux(
        input_names,
        input_names,
        nothing,
        (i::NamedTuple, parameters::Nothing=nothing) -> i
    )
end