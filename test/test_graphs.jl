# using Graphs

# # 创建一个简单的图
# g = SimpleDiGraph(4)
# add_edge!(g, 1, 2)
# add_edge!(g, 2, 3)
# add_edge!(g, 3, 4)

# function upstream_nodes!(upstream_node, graph::SimpleDiGraph, target::Int)
#     in_nodes = inneighbors(graph, target)
#     while !isempty(in_nodes)
#         push!(upstream_node, in_nodes)
#         for n in in_nodes
#             upstream_nodes!(upstream_node, graph, n)
#         end
#     end
#     return upstream_node
# end

# upstream_node = []
# # 获取节点 4 的上游节点
# upstream_nodes!(upstream_node, g, 4)

using Graphs

# 定义一个递归函数来遍历入边邻居
function recursive_inneighbors(graph::SimpleDiGraph, node::Int)
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

# 创建一个图并添加一些边
g = DiGraph(5)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 4)
add_edge!(g, 4, 5)

# 获取节点 3 的所有入边邻居
all_inneighbors = recursive_inneighbors(g, 5)
println(all_inneighbors)  # 输出: Set([1, 2])