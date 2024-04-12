function build_topology(node_names::AbstractVector{Symbol}, reaches::AbstractVector{Pair})
    node_tuple = namedtuple(node_names, collect(1:length(node_names)))
    topology = SimpleDiGraph(length(nodes))
    # add edge
    for reach in reaches
        add_edge!(topology, node_tuple[reach[1]], node_tuple[reach[2]])
    end
    topology
end

#! 自定义了一个找edge index的方法, 在graphs.jl 2.0后更新 https://github.com/JuliaGraphs/Graphs.jl/issues/146
function edge_index(g, src, dst)
    for (idx, edge) in enumerate(edges(g))
        if (edge.src == src) & (edge.dst == dst)
            return idx
        end
    end
    return "not found"
end

function get_all_upstream_node(graph::SimpleDiGraph, node::Any)
    up_node = []
    visited = Set()

    # 辅助递归函数
    function visit(node::Int)
        if node in visited
            return
        end
        push!(visited, node)
        for neighbor in inneighbors(graph, node)
            push!(up_node, neighbor)
            visit(neighbor)  # 递归调用
        end
    end
    # 开始递归
    visit(node)
    return up_node
end