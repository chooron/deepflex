

"""
32 64 128
16     1
8  4   2
"""
d8_codes = [1, 2, 4, 8, 16, 32, 64, 128]
d8_nn_pads = [(1, 1, 2, 0), (2, 0, 2, 0), (2, 0, 1, 1), (2, 0, 0, 2), (1, 1, 0, 2), (0, 2, 0, 2), (0, 2, 1, 1), (0, 2, 2, 0),]


#! nodeinfo: node_id => row_idx, col_idx
function cal_grid_sort(nodeinfo::NamedTuple, flwacc::AbstractMatrix)
    #* 先转换为vec类型
    acc_num = reduce(*, size(flwacc))
    flwacc_vec = vec(reshape(flwacc, 1, acc_num))
    sortperm_vec = sortperm(flwacc_vec)
    #* 按大小计算顺序
    flwacc_sortidx = similar(sortperm_vec)
    for i in 1:length(flwacc_sortidx)
        flwacc_sortidx[sortperm_vec[i]] = i
    end
    #* 生成每个节点的计算顺序
    sortidx_mat = reshape(flwacc_sortidx, size(flwacc)...)
    #* 节点与之对应计算次序的namedtuple
    node_pairs = map(zip(keys(nodeinfo), collect(nodeinfo))) do (k, v)
        k => sortidx_mat[v...]
    end
    #* 根据计算次序对ntp重新排序
    NamedTuple(sort(node_pairs, by=x -> x[2]))
end

function cal_grid_relation(nodeinfo::NamedTuple, flwdir::AbstractMatrix; return_graph=false)
    #* d8的一些基础信息
    d8_direction = [32, 64, 128, 16, 1, 8, 4, 2]
    d8_moveidx = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    d8_dict = Dict([dir => midx for (dir, midx) in zip(d8_direction, d8_moveidx)])

    flwdir_size = size(flwdir)
    di_graph = SimpleDiGraph(length(nodeinfo))
    node_dir_vec = collect(nodeinfo)
    position_dict = Dict([v => i for (i, v) in enumerate(node_dir_vec)])
    #* 根据flwdir构建graphs
    for row_idx in 1:size(flwdir)[1]
        for col_idx in 1:size(flwdir)[2]
            if !(flwdir[row_idx, col_idx] in d8_direction)
                continue
            end
            move_idx_ = d8_dict[flwdir[row_idx, col_idx]]
            target_row_idx, target_col_idx = row_idx + move_idx_[1], col_idx + move_idx_[2]
            if (target_row_idx < 1) | (target_row_idx > flwdir_size[1])
                continue
            end
            if (target_col_idx < 1) | (target_col_idx > flwdir_size[2])
                continue
            end
            add_edge!(di_graph, position_dict[(row_idx, col_idx)], position_dict[(target_row_idx, target_col_idx)])
        end
    end

    if return_graph
        di_graph
    else
        #* 获取graph所有点
        all_vertices = Graphs.vertices(di_graph)
        #* 获取每个点的matidx和其对应的输入点的mat idx
        inflow_relation = NamedTuple(map(all_vertices) do vertice
            keys(nodeinfo)[vertice] => [keys(nodeinfo)[up_vert] for up_vert in Graphs.inneighbors(di_graph, vertice)]
        end)
        inflow_relation
    end
end

function build_grid_graph(nodeinfo::NamedTuple, flwdir::AbstractMatrix)
    d8_direction = [32, 64, 128, 16, 1, 8, 4, 2]
    d8_moveidx = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    d8_dict = Dict([dir => midx for (dir, midx) in zip(d8_direction, d8_moveidx)])

    flwdir_size = size(flwdir)
    di_graph = SimpleDiGraph(length(nodeinfo))
    node_dir_vec = collect(nodeinfo)
    position_dict = Dict([v => i for (i, v) in enumerate(node_dir_vec)])
    #* 根据flwdir构建graphs
    for row_idx in 1:size(flwdir)[1]
        for col_idx in 1:size(flwdir)[2]
            if !(flwdir[row_idx, col_idx] in d8_direction)
                continue
            end
            move_idx_ = d8_dict[flwdir[row_idx, col_idx]]
            target_row_idx, target_col_idx = row_idx + move_idx_[1], col_idx + move_idx_[2]
            if (target_row_idx < 1) | (target_row_idx > flwdir_size[1])
                continue
            end
            if (target_col_idx < 1) | (target_col_idx > flwdir_size[2])
                continue
            end
            add_edge!(di_graph, position_dict[(row_idx, col_idx)], position_dict[(target_row_idx, target_col_idx)])
        end
    end
    di_graph
end

function get_route_len(graph::Graphs.DiGraph, target::Int)
    route_len_dict = Dict()
    for vt in Graphs.vertices(di_graph)
        route_len_dict[vt] = length(Graphs.a_star(graph, vt, target))
    end
    route_len_dict
end


function grid_routing(flow::AbstractMatrix, flwdir::AbstractMatrix)
    flow_routed = sum(collect([pad_zeros(flow .* (flwdir .== code), arg)
                               for (code, arg) in zip(d8_codes, d8_nn_pads)]))

    flow_routed[2:size(flow)[1]+1, 2:size(flow)[2]+1]
end

function grid_routingv2(flow::AbstractMatrix, flwdir::AbstractMatrix)
    flow_routed = sum(collect([pad(flow .* (flwdir .== code), eltype(flow)(0.0), arg1, arg2)
                               for (code, arg1, arg2) in zip(d8_codes, d8_pads_args1, d8_pads_args2)]))

    flow_routed[2:size(flow)[1]+1, 2:size(flow)[2]+1]
end

function matrix_idx(arr_size::Tuple, idx::Tuple)
    n = arr_size[2]
    (idx[1] - 1) * n + idx[2]
end
