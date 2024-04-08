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