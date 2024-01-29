@kwdef struct Flow{T<:Number} <: AbstractFunc
    input_names::Vector{Symbol}
    output_names::Vector{Symbol} = [:Flow]
    parameters::ComponentVector{T}
    weight::ComponentVector{T}
end

function Flow(input_names::Vector{Symbol}; parameters::ComponentVector{T}, weights::ComponentVector{T}) where {T<:Number}
    Surfaceflow{T}(input_names=input_names, parameters=parameters, weights=weights)
end

function get_output(ele::Flow; input::ComponentVector{T}) where {T<:Number}
    args = [input[input_nm] for input_nm in ele.input_names]
    ComponentVector(; Dict(first(ele.output_names) => surfaceflow.(args...; ele.parameters...))...)
end

function flow(baseflow::T, surfaceflow::T) where {T<:Number}
    baseflow + surfaceflow
end