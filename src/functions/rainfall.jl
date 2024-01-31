@kwdef struct Rainfall{T<:Number} <: AbstractFunc
    input_names::Vector{Symbol}
    output_names::Vector{Symbol} = [:Rainfall]
    parameters::ComponentVector{T}
end

function Rainfall(input_names::Vector{Symbol}; parameters::ComponentVector{T}) where {T<:Number}
    Rainfall{T}(input_names=input_names, parameters=parameters)
end

function get_output(ele::Rainfall; input::ComponentVector{T}) where {T<:Number}
    args = [input[input_nm] for input_nm in ele.input_names]
    ComponentVector(; Dict(first(ele.output_names) => rainfall.(args...; ele.parameters...))...)
end

function rainfall(Prcp::T, Temp::T; Tmin::T) where {T<:Number}
    step_func(Temp - Tmin) * Prcp
end
