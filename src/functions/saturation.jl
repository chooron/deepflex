@kwdef struct Saturation{T<:Number} <: AbstractFunc
    input_names::Vector{Symbol}
    output_names::Vector{Symbol} = [:Rainfall]
    parameters::ComponentVector{T}
end

function Saturation(input_names::Vector{Symbol}; parameters::ComponentVector{T}) where {T<:Number}
    Saturation{T}(input_names=input_names, parameters=parameters)
end

function get_output(ele::Saturation; input::ComponentVector{T}) where {T<:Number}
    args = [input[input_nm] for input_nm in ele.input_names]
    ComponentVector(; Dict(first(ele.output_names) => rainfall.(args...; ele.parameters...))...)
end

function saturation(SoilWater::T, Rainfall::T; x1::T) where {T<:Number}
    Rainfall * (1 - (SoilWater / x1)^2)
end
