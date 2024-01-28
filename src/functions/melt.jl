@kwdef struct Melt{T<:Number} <: AbstractFunc
    input_names::Vector{Symbol}
    output_name::Symbol = :Melt
    parameters::Dict{Symbol,T}
end

function Melt(input_names::Vector{Symbol}; parameters::Dict{Symbol,T}) where {T<:Number}
    Melt{T}(input_names=input_names, parameters=parameters)
end

function get_output(ele::Melt; input::ComponentVector{T}) where {T<:Number}
    args = [input[input_nm] for input_nm in ele.input_names]
    ComponentVector(; Dict(ele.output_name => melt.(args...; ele.parameters...))...)
end

function melt(SoilWater::T, Temp::T; Tmax::T, Df::T) where {T<:Number}
    step_func(Temp - Tmax) * step_func(SoilWater) * min(SoilWater, Df * (Temp - Tmax))
end