@kwdef struct Flow{T<:Number} <: AbstractFunc
    input_names::Vector{Symbol}
    output_names::Vector{Symbol} = [:Flow]
    parameters::ComponentVector{T}
end

function Flow(input_names::Vector{Symbol}; parameters::ComponentVector{T}) where {T<:Number}
    Flow{T}(input_names=input_names, parameters=parameters)
end

function get_output(ele::Flow; input::ComponentVector{T}) where {T<:Number}
    args = [input[input_nm] for input_nm in ele.input_names]
    ComponentVector(; Dict(first(ele.output_names) => flow.(args...; ele.parameters...))...)
end

function flow(Baseflow::T, Surfaceflow::T) where {T<:Number}
    Baseflow + Surfaceflow
end