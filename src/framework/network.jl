mutable struct Network{N<:Node,T<:Number} <: Component
    id::String
    topology::AbstractGraph
end

function Network(; id::String, nodes::Dict{Symbol,N}, topology::AbstractGraph)
    topology = MetaDiGraph(topology)
    for (name, node) in nodes
        set_props!(topology, name, Dict(:node => node))
    end
    Network(id, topology)
end

function get_output(network::N, input::Dict{Symbol,Dict{Symbol,Vector{T}}}) where {N<:Network,T<:Number}
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
        tmp_ouput = get_output(tmp_node, input=input[node_nm])
        for up_node_nm in inneighbors(network.topology, node_nm)
            tmp_up_node = get_prop(network.topology, up_node_nm, :node)
            # calculate upstream subbasin routing output
            routing_out = tmp_up_node.routing_func(tmp_up_node.fluxes, tmp_up_node.target_names) * total_area[up_node_nm]
            # combine output and routing out
            for (key, value) in tmp_ouput
                tmp_ouput[key] = value * tmp_node.area + routing_out[key] * node_nm.area
            end
        end
        # Divide the total output result by the total area
        for (key, value) in tmp_ouput
            tmp_ouput[key] /= node_nm.total_area
        end
        output[node_nm] = tmp_ouput
    end
    return output
end