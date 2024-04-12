using Graphs
using MetaGraphs
# using GraphPlot

# 创建一个空图
g = DiGraph(5)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 2, 4)
add_edge!(g, 3, 5)
add_edge!(g, 4, 5)

inneighbors(g, 4)
inneighbors(g, 3)
outneighbors(g, 3)

Graphs.Edge(2,4)

# 添加不同类型的节点
Graphs.add_vertices!(g, "A")
add_vertex!(g, 2)
add_vertex!(g, 3.5, label="FloatNode")

# 添加带有类型的边
add_edge!(g, "A", 2, label="SymbolToIntegerEdge")
add_edge!(g, 2, 4)

# 可视化图形
# gplot(g, nodelabel=1:nv(g), edgelabel=1:ne(g), nodefillc="lightblue", nodeshape=:circle, nodefontsize=10, edgetype=:arrow)
# 提取出每个