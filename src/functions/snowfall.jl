@kwdef struct Snowfall{T<:Number} <: AbstractFunc
    input_names::Vector{Symbol}
    output_name::Symbol = :Snowfall
    parameters::Dict{Symbol,T}
end

function Snowfall(input_names::Vector{Symbol}; parameters::Dict{Symbol,T}) where {T<:Number}
    Snowfall{T}(id=id, input_names=input_names, parameters=parameters)
end

function get_output(ele::Snowfall; input::ComponentVector{T}) where {T<:Number}
    args = [input[input_nm] for input_nm in ele.input_names]
    ComponentVector(; Dict(ele.output_name => snowfall.(args...; ele.parameters...))...)
end

function snowfall(Pcrp::T, Temp::T, Tmin::T) where {T<:Number}
    step_func(Tmin - Temp) * Pcrp
end