struct RiverNetwork <: AbstractRiverNetwork
    name::Symbol
    #* 水文计算节点,出流量为discharge(m^3/s)
    nodes::AbstractVector{HydroNode}
    #* 河道对象,用于洪水演进模拟,包含基于水文的马斯京根模拟和一维水动计算方式(todo) 
    reaches::AbstractVector{AbstractReach}
    #* 河网
    topology::AbstractGraph
end

struct GridNetwork <: AbstractGridNetwork
    # todo 通过网格实现flux传播
end

function RiverNetwork(;
    name::Symbol,
    nodes::AbstractVector{HydroNode},
    reaches::AbstractVector{AbstractReach}
)
    node_names = [node.name for node in name]
    reach_nodes = [reach.upstream => reach.downstream for reach in reaches]
    topology = build_topology(node_names, reach_nodes)
    RiverNetwork(name, nodes, reaches, topology)
end

function (network::RiverNetwork)(
    input::NamedTuple,
    params::NamedTuple,
    init_states::NamedTuple;
    dt::Number=1.0,
)
    # 先同时计算各个节点，然后自上而下计算河网汇流
    # 1.计算然后收集计算结果
    # 2.按层次遍历所有节点，判断这个节点上游最近一层的所有节点，找到该节点与上游节点的所有边，根据边获取所有的演进公式并带入计算
    node_names = [node.nm for nm in network.nodes]
    # todo 添加node的并行计算
    node_result_tuple = namedtuple(
        node_names, [node(input, params, init_states, step=false)[:discharge] for nm in network.nodes]
    )

    for node_idx in topological_sort(network.topology)
        node_name = node_names[node_idx]
        if length(get_all_upstream_node(network.topology, node_nm)) != 0
            node_result = node_result_tuple[node_name]
            for up_node_idx in get_all_upstream_node(network.topology, node_nm)
                up_node_name = node_names[up_node_idx]
                up_node_result = node_result_tuple[up_node_name]
                node_result = node_result .+ network.reach[edge_index(network.topology, up_node_idx, idx)](up_node_result, dt)
            end
        end
        node_result_tuple = merge(node_result_tuple, namedtuple([node_name], [node_result_tuple[node_name]]))
    end
    node_result_tuple
end


function (network::GridNetwork)(
    input::NamedTuple,
    params::NamedTuple,
    init_states::NamedTuple
)
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