using Graphs
using MetaGraphsNext

colors = MetaGraph(
    Graph();  # underlying graph structure
    label_type=Symbol,  # color name
    vertex_data_type=NTuple{3,Int},  # RGB code
    edge_data_type=Symbol,  # result of the addition between two colors
    graph_data="additive colors",  # tag for the whole graph
)

colors[:red] = (255, 0, 0);
colors[:green] = (0, 255, 0);
colors[:blue] = (0, 0, 255);

colors[:red, :green] = :yellow;
colors[:red, :blue] = :magenta;
colors[:green, :blue] = :cyan;

graph = Graph(Edge.([(1, 2), (1, 3), (2, 3)]))
vertices_description = [:red => (255, 0, 0), :green => (0, 255, 0), :blue => (0, 0, 255)]
edges_description = [
    (:red, :green) => :yellow, (:red, :blue) => :magenta, (:green, :blue) => :cyan
]

colors2 = MetaGraph(graph, vertices_description, edges_description, "additive colors")
colors2 == colors