@kwdef struct Evap{T<:Number} <: AbstractFunc
    input_names::Vector{Symbol}
    output_names::Vector{Symbol} = [:Evap]
    parameters::ComponentVector{T}
end

function Evap(input_names::Vector{Symbol}; parameters::ComponentVector{T}) where {T<:Number}
    Evap{T}(input_names=input_names, parameters=parameters)
end

function get_output(ele::Evap; input::ComponentVector{T}) where {T<:Number}
    args = [input[input_nm] for input_nm in ele.input_names]
    ComponentVector(; Dict(first(ele.output_names) => evap.(args...; ele.parameters...))...)
end

function evap(SoilWater::T, Pet::T; Smax::T) where {T<:Number}
    step_func(SoilWater) * step_func(SoilWater - Smax) * Pet + step_func(SoilWater) * step_func(Smax - SoilWater) * Pet * (SoilWater / Smax)
end

function evap(SoilWater::T, Prcp::T, Pet::T; x1::T) where {T<:Number}
    step_func(Pet - Prcp) * (Pet - Prcp) * (2 * SoilWater / x1 - (SoilWater / x1)^2)
end