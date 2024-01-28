@kwdef struct Rainfall{T<:Number} <: AbstractFunc
    input_names::Vector{Symbol}
    output_name::Symbol = :Rainfall
    parameters::Dict{Symbol,T}
end

function Rainfall(input_names::Vector{Symbol}; parameters::Dict{Symbol,T}) where {T<:Number}
    Rainfall{T}(input_names=input_names, parameters=parameters)
end

function get_output(ele::Rainfall; input::ComponentVector{T}) where {T<:Number}
    args = [input[input_nm] for input_nm in ele.input_names]
    ComponentVector(;Dict(ele.output_name=>rainfall.(args...; ele.parameters...))...)
end

function rainfall(Prcp::T, Temp::T; Tmin::T) where{T<:Number}
    step_func(Temp - Tmin) * Prcp
end
