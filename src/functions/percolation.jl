@kwdef struct Percolation{T<:Number} <: AbstractFunc
    input_names::Vector{Symbol}
    output_names::Vector{Symbol} = [:Melt]
    parameters::ComponentVector{T}
end

function Percolation(input_names::Vector{Symbol}; parameters::ComponentVector{T}) where {T<:Number}
    Percolation{T}(input_names=input_names, parameters=parameters)
end

function get_output(ele::Percolation; input::ComponentVector{T}) where {T<:Number}
    args = [input[input_nm] for input_nm in ele.input_names]
    ComponentVector(; Dict(first(ele.output_names) => melt.(args...; ele.parameters...))...)
end

function percolation(SoilWater::T; x1::T) where {T<:Number}
    (x1^(-4)) / 4 * ((4 / 9)^(-4)) * (SoilWater^5)
end