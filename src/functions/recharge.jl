@kwdef struct Recharge{T<:Number} <: AbstractFunc
    input_names::Vector{Symbol}
    output_names::Vector{Symbol} = [:Rainfall]
    parameters::ComponentVector{T}
end

function Recharge(input_names::Vector{Symbol}; parameters::ComponentVector{T}) where {T<:Number}
    Recharge{T}(input_names=input_names, parameters=parameters)
end

function get_output(ele::Recharge; input::ComponentVector{T}) where {T<:Number}
    args = [input[input_nm] for input_nm in ele.input_names]
    ComponentVector(; Dict(first(ele.output_names) => rainfall.(args...; ele.parameters...))...)
end

function recharge(RoutingStore::T; x2::T, x3::T, ω::T) where {T<:Number}
    x2 / (x3^ω) * RoutingStore^ω
end