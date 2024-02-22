function Splitter(input_names::Vector{Symbol}; parameters::ComponentVector{T}) where {T<:Number}
    SimpleFlux(
        input_names,
        collect(keys(parameters)),
        parameters,
        splitter_func
    )
end

function splitter_func(
    input::NamedTuple,
    parameters::NamedTuple
)::NamedTuple
    tmp_input = input(first(input_names))
    (;Dict(k=parameters[k] .* tmp_input for k in keys(parameters))...)
end
