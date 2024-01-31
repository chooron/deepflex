@kwdef struct Snowfall{T<:Number} <: AbstractFunc
    input_names::Vector{Symbol}
    output_names::Vector{Symbol} = [:Snowfall]
    parameters::ComponentVector{T}
end

function Snowfall(input_names::Vector{Symbol}; parameters::ComponentVector{T}) where {T<:Number}
    Snowfall{T}(input_names=input_names, parameters=parameters)
end

function get_output(ele::Snowfall; input::ComponentVector{T}) where {T<:Number}
    args = [input[input_nm] for input_nm in ele.input_names]
    ComponentVector(; Dict(first(ele.output_names) => snowfall.(args...; ele.parameters...))...)
end

function snowfall(Pcrp::T, Temp::T; Tmin::T) where {T<:Number}
    step_func(Tmin - Temp) * Pcrp
end