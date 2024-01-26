@kwdef struct Melt{T<:Number} <: ParameterizedElement
    id::String
    input_names::Vector{Symbol}
    output_names::Vector{Symbol} = [:Qb]
    parameters::Dict{Symbol,T}
end

function Melt(; id::String, input_names::Vector{Symbol}, parameters::Dict{Symbol,T}) where {T<:Number}
    Melt{T}(id=id, input_names=input_names, parameters=parameters)
end

function get_output(ele::Melt; input::ComponentVector{T}) where {T<:Number}
    args = [input[input_nm] for input_nm in ele.input_names]
    melt.(args...; ele.parameters...)
end

function melt(S::T, Temp::T; Tmax::T, Df::T) where {T<:Number}
    step_func(Temp - Tmax) * step_func(S) * min(S, Df * (Temp - Tmax))
end