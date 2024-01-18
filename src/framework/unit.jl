function get_output(unit::Unit; input::Dict{Symbol,Vector{T}}) where {T<:Number}
    # initialize unit fluxes
    unit.fluxes = input
    # traversal of the directed graph
    for idx in topological_sort(unit.topology)
        tmp_ele = get_prop(unit.topology, idx, :ele)
        tmp_fluxes = get_output(tmp_ele, input=unit.fluxes)
        merge!(unit.fluxes, tmp_fluxes)
    end
    return unit.fluxes
end

