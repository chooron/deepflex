@kwdef struct Snowfall{T<:Number} <: ParameterizedElement
    id::String
    input_names::Vector{Symbol}
    output_names::Vector{Symbol} = [:Pet]
    parameters::Dict{Symbol,T}
end

function Snowfall(; id::String, input_names::Vector{Symbol}, parameters::Dict{Symbol,T}) where {T<:Number}
    Snowfall{T}(id=id, input_names=input_names, parameters=parameters)
end

function get_output(ele::Snowfall; input::ComponentVector{T}) where {T<:Number}
    args = [input[input_nm] for input_nm in ele.input_names]
    snowfall.(args...; ele.parameters...)
end

function snowfall(Pcrp::T, Temp::T, Tmin::T) where {T<:Number}
    step_func(Tmin - Temp) * Pcrp
end