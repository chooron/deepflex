function Splitter(input_names::Vector{Symbol}; parameters::ComponentVector{T}) where {T<:Number}
    SimpleFlux{T}(
        input_names,
        collect(keys(parameters)),
        parameters,
        splitter_func
    )
end

function splitter_func(
    input::ComponentVector{T},
    parameters::ComponentVector{T}
) where {T<:Number}
    tmp_input = input(first(input_names))
    ComponentVector(;Dict(k=parameters[k] * tmp_input for k in keys(parameters))...)
end
