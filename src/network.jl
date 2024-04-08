struct Network{T<:Number} <: AbstractComponent
    name::Symbol
    nodes::AbstractVector{HydroNode}
    areas::NamedTuple
    routes::NamedTuple
    topology::AbstractGraph
end

function Network(;
    name::Symbol,
    nodes::AbstractVector{HydroNode},
    areas::NamedTuple,
    routes::NamedTuple,
    topology::AbstractGraph
)
    topology = MetaDiGraph(topology)
    for node in nodes
        node_name = node.nameinfo.name
        set_props!(topology, node_name, Dict(:node => node, :area => areas[node_name]))
    end
    # todo 添加连接线的属性
    Network(name, nodes, topology)
end

function (network::N)(
    input::NamedTuple,
    params::NamedTuple,
    init_states::NamedTuple
) where {N<:Network,T<:Number}
    #! 这个图貌似不需要添加属性进去
    # calculate subbasin and it's all upstream total area
    total_area::Dict{Symbol,T} = Dict()
    output::Dict{Symbol,Dict{Symbol,Vector{T}}} = Dict()

    for node_nm in topological_sort(network.topology)
        tmp_area = get_prop(network.topology, node_nm, :node).area
        for up_node_nm in get_all_upstream_node(network.topology, node_nm)
            tmp_area += get_prop(network.topology, up_node_nm, :node).area
        end
        total_area[node_nm] = tmp_area
    end

    for node_nm in topological_sort(network.topology)
        tmp_node = get_prop(network.topology, node_nm, :node)
        # calculate current subbasin output
        tmp_ouput = tmp_node(input, params, init_states)
        for up_node_nm in inneighbors(network.topology, node_nm)
            tmp_up_node = get_prop(network.topology, up_node_nm, :node)
            # calculate upstream subbasin routing output
            routing_out = tmp_up_node.routing_func(tmp_up_node.fluxes, tmp_up_node.target_names) * total_area[up_node_nm]
            # combine output and routing out
            for (key, value) in pairs(tmp_ouput)
                tmp_ouput[key] = value * tmp_node.area + routing_out[key] * node_nm.area
            end
        end
        # Divide the total output result by the total area
        for (key, value) in pairs(tmp_ouput)
            tmp_ouput[key] /= node_nm.total_area
        end
        output[node_nm] = tmp_ouput
    end
    return output
end