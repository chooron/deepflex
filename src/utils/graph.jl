function build_compute_topology(fluxes::AbstractVector{<:AbstractFlux})
    # 构建函数之间的计算图
    input_names, output_names = get_func_io_names(fluxes)
    input_names_ntp = namedtuple(input_names, collect(1:length(input_names)))
    output_names_ntp = namedtuple(output_names, collect(1:length(output_names)))

    topology = SimpleDiGraph(length(nodes))
    for flux in fluxes
        tmp_input_names, tmp_output_names = get_input_names(flux), get_output_names(flux)
        for ipnm in tmp_input_names
            for opnm in tmp_output_names
                add_edge!(topology, input_names_ntp[ipnm], output_names_ntp[opnm])
            end
        end
    end
    topology
end


function build_river_topology(node_names::AbstractVector{Symbol}, reaches::AbstractVector)
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

