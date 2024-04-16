struct RiverNetwork <: AbstractRiverNetwork
    name::Symbol
    input_names::NamedTuple
    output_names::NamedTuple
    state_names::NamedTuple
    param_names::NamedTuple
    #* 水文计算节点,出流量为discharge(m^3/s)
    nodes::AbstractVector{HydroNode}
    #* 河道对象,用于洪水演进模拟,包含基于水文的马斯京根模拟和一维水动计算方式(todo) 
    reaches::AbstractVector{AbstractReach}
    #* 河网
    topology::AbstractGraph
end

struct GridNetwork <: AbstractGridNetwork
    # todo 通过网格实现flux传播
    # 如何根据网格的高程分布计算出网络的汇流关系
end

function RiverNetwork(;
    name::Symbol,
    nodes::AbstractVector{HydroNode},
    reaches::AbstractVector{AbstractReach}
)
    input_names, output_names, state_names, param_names = get_network_infos(nodes, reaches)
    node_names = [node.name for node in nodes]
    reach_nodes = [reach.upstream => reach.downstream for reach in reaches]
    topology = build_topology(node_names, reach_nodes)
    RiverNetwork(
        name,
        input_names,
        output_names,
        state_names,
        param_names,
        nodes,
        reaches,
        topology
    )
end

function get_network_infos(nodes::AbstractVector{HydroNode}, reaches::AbstractVector{AbstractReach})
    input_names, output_names, state_names, param_names = NamedTuple(), NamedTuple(), NamedTuple(), NamedTuple()
    #* 获取nodes信息
    for node in network.nodes
        input_names = merge(input_names, namedtuple([node.name], [node.input_names]))
        output_names = merge(output_names, namedtuple([node.name], [node.output_names]))
        state_names = merge(state_names, namedtuple([node.name], [node.state_names]))
        param_names = merge(param_names, namedtuple([node.name], [node.param_names]))
    end
    #* 获取reach的参数信息
    for reach in network.reaches
        #! 暂时以马斯京根作为河道的主流计算方式
        param_names = merge(param_names, namedtuple([reach.name], [:k, :x]))
    end
    input_names, output_names, state_names, param_names
end

function (network::RiverNetwork)(
    input::NamedTuple,
    params::ComponentVector,
    init_states::ComponentVector;
    dt::Number=1.0,
)
    # 先同时计算各个节点，然后自上而下计算河网汇流
    # 1.计算然后收集计算结果
    # 2.按层次遍历所有节点，判断这个节点上游最近一层的所有节点，找到该节点与上游节点的所有边，根据边获取所有的演进公式并带入计算
    node_names = [node.name for node in network.nodes]
    # todo 添加node的并行计算
    node_result_tuple = namedtuple(
        node_names, [node(input, params, init_states, step=false)[:discharge] for nm in network.nodes]
    )
    for node_idx in topological_sort(network.topology)
        node_name = node_names[node_idx]
        node_result = node_result_tuple[node_name]
        for up_node_idx in inneighbors(network.topology, node_idx)
            up_node_name = node_names[up_node_idx]
            up_node_result = node_result_tuple[up_node_name]
            node_result = node_result .+ network.reaches[edge_index(network.topology, up_node_idx, node_idx)](up_node_result, dt)
        end
        node_result_tuple = merge(node_result_tuple, namedtuple([node_name], [node_result]))
    end
    node_result_tuple
end

function (network::GridNetwork)(
    input::NamedTuple,
    params::ComponentVector,
    init_states::ComponentVector
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