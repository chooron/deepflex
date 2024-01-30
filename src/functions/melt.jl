@kwdef struct Melt{T<:Number} <: AbstractFunc
    input_names::Vector{Symbol}
    output_names::Vector{Symbol} = [:Melt]
    parameters::ComponentVector{T}
    weights::ComponentVector{T}
end

function Melt(input_names::Vector{Symbol}; parameters::ComponentVector{T}, weights::ComponentVector{T}) where {T<:Number}
    Melt{T}(input_names=input_names, parameters=parameters, weights=weights)
end

function get_output(ele::Melt; input::ComponentVector{T}) where {T<:Number}
    args = [input[input_nm] for input_nm in ele.input_names]
    ComponentVector(; Dict(first(ele.output_names) => melt.(args...; ele.parameters...))...)
end

function melt(SnowWater::T, Temp::T; Tmax::T, Df::T) where {T<:Number}
    step_func(Temp - Tmax) * step_func(SnowWater) * min(SnowWater, Df * (Temp - Tmax))
end