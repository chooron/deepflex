@kwdef struct Baseflow{T<:Number} <: AbstractFunc
    input_names::Vector{Symbol}
    output_names::Vector{Symbol} = [:bf]
    parameters::ComponentVector{T}
end

function Baseflow(input_names::Vector{Symbol}; parameters::ComponentVector{T}) where {T<:Number}
    Baseflow{T}(input_names=input_names, parameters=parameters)
end

function (::Baseflow)(input::ComponentVector{T}) where {T<:Number}
    ComponentVector(; Dict(first(ele.output_names) => baseflow.(input; ele.parameters))...)
end

function get_output(ele::Baseflow; input::ComponentVector{T}) where {T<:Number}
    args = [input[input_nm] for input_nm in ele.input_names]
    ComponentVector(; Dict(first(ele.output_names) => baseflow.(args...; ele.parameters...))...)
end

function baseflow(
    input::ComponentVector{T,Vector{T},Tuple{Axis{(sw = 1)}}},
    params::ComponentVector{T,Vector{T},Tuple{Axis{(Smax=1, Qmax=2, f=3)}}}
) where {T<:Number}
    @variables t sw(t) bf(t)
    @parameters Smax Qmax f
    bf ~ step_func(sw) * step_func(sw - Smax) * Qmax + step_func(sw) * step_func(Smax - sw) * Qmax * exp(-f * (Smax - sw))
end

function baseflow(SoilWater::T; Smax::T, Qmax::T, f::T) where {T<:Number}
    step_func(SoilWater) * step_func(SoilWater - Smax) * Qmax + step_func(SoilWater) * step_func(Smax - SoilWater) * Qmax * exp(-f * (Smax - SoilWater))
end

function baseflow(input::ComponentVector{T,Vector{T},Tuple{Axis{(sw = 1)}}},
    params::ComponentVector{T,Vector{T},Tuple{Axis{(Smax=1, Qmax=2, f=3)}}}
) where {T<:Number}
    @variables t rs bf
    @parameters x3 γ
    x3^(1 - γ) / (γ - 1) * rs^γ
end
