@kwdef struct Splitter{T<:Number} <: AbstractFunc
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
    parameters::ComponentVector{T}
end

function Splitter(input_names::Vector{Symbol}; parameters::ComponentVector{T}) where {T<:Number}
    Splitter{T}(input_names=input_names, output_names=collect(keys(parameters)), parameters=parameters)
end

function get_output(ele::Splitter; input::ComponentVector{T}) where {T<:Number}
    tmp_input = input(first(input_names))
    ComponentVector(; Dict(k => ele.parameters[k] * tmp_input)...)
end
