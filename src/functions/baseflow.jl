@kwdef struct Baseflow{T<:Number} <: AbstractFunc
    input_names::Vector{Symbol}
    output_names::Vector{Symbol} = [:Baseflow]
    parameters::ComponentVector{T}
    weights::ComponentVector{T}
end

function Baseflow(input_names::Vector{Symbol}; parameters::ComponentVector{T}, weights::ComponentVector{T}) where {T<:Number}
    Baseflow{T}(input_names=input_names, parameters=parameters, weights=weights)
end

function get_output(ele::Baseflow; input::ComponentVector{T}) where {T<:Number}
    args = [input[input_nm] for input_nm in ele.input_names]
    ComponentVector(; Dict(first(ele.output_names) => baseflow.(args...; ele.parameters...))...)
end

function baseflow(SoilWater::T; Smax::T, Qmax::T, f::T) where {T<:Number}
    step_func(SoilWater) * step_func(SoilWater - Smax) * Qmax + step_func(SoilWater) * step_func(Smax - SoilWater) * Qmax * exp(-f * (Smax - SoilWater))
end