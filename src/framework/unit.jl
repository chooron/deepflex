@kwdef mutable struct Unit{T} <: Component where {T<:Number}
    id::String

    # model structure
    stucture::AbstractGraph

    # inner variables
    fluxes::Dict{Symbol,Vector{T}} = Dict()

    # attribute
    input_names::Vector{Symbol} = []
end


function get_output(unit::Unit; input::Dict{Symbol,Vector{T}}) where {T<:Number}
    # initialize unit fluxes
    unit.fluxes = input
    # traversal of the directed graph
    for idx in topological_sort(unit.stucture)
        tmp_ele = get_prop(unit.stucture, idx, :ele)
        tmp_fluxes = get_output(tmp_ele, input=unit.fluxes)
        merge!(unit.fluxes, tmp_fluxes)
    end
    return unit.fluxes
end

