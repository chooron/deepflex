using Graphs
using MetaGraphsNext
using MetaGraphs

# colors = MetaGraph(
#     DiGraph();  # underlying graph structure
#     label_type=Symbol,  # color name
#     vertex_data_type=NTuple{3,Int},  # RGB code
#     edge_data_type=Symbol,  # result of the addition between two colors
#     graph_data="additive colors",  # tag for the whole graph
# )

# colors[:red] = (255, 0, 0);
# colors[:green] = (0, 255, 0);
# colors[:blue] = (0, 0, 255);

# colors[:red, :green] = :yellow;
# colors[:red, :blue] = :magenta;
# colors[:green, :blue] = :cyan;


# for g in topological_sort(colors)
#     println(MetaGraphsNext.(colors, 2, :rgb_code))
# end

# collect(labels(colors))


dag = SimpleDiGraph(4)
add_edge!(dag, 1, 2)
add_edge!(dag, 2, 3)
add_edge!(dag, 3, 4)
topology = MetaDiGraph(SimpleDiGraph(4))
for idx in topological_sort(topology)
    println(idx)
    # tmp_fluxes = get_output(tmp_ele, unit.fluxes)
    # merge!(unit.fluxes, tmp_fluxes)
end
for idx in topological_sort(dag)
    println(idx)
    # tmp_fluxes = get_output(tmp_ele, unit.fluxes)
    # merge!(unit.fluxes, tmp_fluxes)
end
# colors_unstable = MetaGraphs.MetaGraph(DiGraph(3))
# MetaGraphs.set_indexing_prop!(colors_unstable, :label)

# MetaGraphs.set_prop!(colors_unstable, :graph_tag, "additive colors")

# MetaGraphs.set_props!(colors_unstable, 1, Dict(:label => :red, :rgb_code => (255, 0, 0)))
# MetaGraphs.set_props!(colors_unstable, 2, Dict(:label => :green, :rgb_code => (0, 255, 0)))
# MetaGraphs.set_props!(colors_unstable, 3, Dict(:label => :blue, :rgb_code => (0, 0, 255)))

# MetaGraphs.set_prop!(colors_unstable, 1, 2, :addition_result, :yellow)
# MetaGraphs.set_prop!(colors_unstable, 1, 3, :addition_result, :magenta)
# MetaGraphs.set_prop!(colors_unstable, 2, 3, :addition_result, :cyan);

# for g in topological_sort(colors)
#     println(MetaGraphs.get_prop(colors_unstable, g, :rgb_code))
# end
