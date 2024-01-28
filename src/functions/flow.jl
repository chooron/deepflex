@kwdef struct Flow{T<:Number} <: AbstractFunc
    input_names::Vector{Symbol}
    output_name::Symbol = :Flow
    parameters::Dict{Symbol,T}
end

function Flow(input_names::Vector{Symbol}; parameters::Dict{Symbol,T}) where {T<:Number}
    Surfaceflow{T}(input_names=input_names, parameters=parameters)
end

function get_output(ele::Flow; input::ComponentVector{T}) where {T<:Number}
    args = [input[input_nm] for input_nm in ele.input_names]
    ComponentVector(; Dict(ele.output_name => surfaceflow.(args...; ele.parameters...))...)
end

function flow(baseflow::T, surfaceflow::T) where {T<:Number}
    baseflow + surfaceflow
end