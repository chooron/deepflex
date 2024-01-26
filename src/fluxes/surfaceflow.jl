@kwdef struct Surfaceflow{T<:Number} <: ParameterizedElement
    input_names::Vector{Symbol}
    output_names::Vector{Symbol} = [:Pet]
    parameters::Dict{Symbol,T}
end

function Surfaceflow(; input_names::Vector{Symbol}, parameters::Dict{Symbol,T}) where {T<:Number}
    Surfaceflow{T}(id=id, input_names=input_names, parameters=parameters)
end

function get_output(ele::Surfaceflow; input::ComponentVector{T}) where {T<:Number}
    args = [input[input_nm] for input_nm in ele.input_names]
    surfaceflow.(args...; ele.parameters...)
end

function surfaceflow(SoilWater::T; Smax::T) where {T<:Number}
    step_func(SoilWater) * step_func(SoilWater .- Smax) .* (SoilWater .- Smax)
end