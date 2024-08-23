

"""
32 64 128
16     1
8  4   2
"""
d8_codes = [1, 2, 4, 8, 16, 32, 64, 128]
d8_nn_pads = [(1, 1, 2, 0), (1, 1, 2, 0), (2, 0, 2, 0), (2, 0, 1, 1), (2, 0, 0, 2), (1, 1, 0, 2), (0, 2, 0, 2), (0, 2, 1, 1), (0, 2, 2, 0),]

"""
$(SIGNATURES)

Construct calculation graphs based on all common hydrological fluxes in hydrological components
"""
function sort_fluxes(fluxes::AbstractVector{<:AbstractComponent})
    input_names = reduce(union, get_input_names.(fluxes))
    output_names = reduce(union, get_output_names.(fluxes))
    input_names = setdiff(input_names, output_names)
    output_names = setdiff(output_names, input_names)

    #* 构建flux输出名称与实例的namedtuple
    fluxes_ntp = reduce(merge, map(fluxes) do flux
        tmp_output_names = get_output_names(flux)
        NamedTuple{Tuple(tmp_output_names)}(repeat([flux], length(tmp_output_names)))
    end)

    #* 构建flux的有向计算图
    var_names = vcat(input_names, output_names)
    var_names_ntp = NamedTuple{Tuple(var_names)}(1:length(var_names))
    digraph = SimpleDiGraph(length(var_names))
    for flux in fluxes
        tmp_input_names, tmp_output_names = get_input_names(flux), get_output_names(flux)
        for ipnm in tmp_input_names
            for opnm in tmp_output_names
                println((ipnm => opnm))
                add_edge!(digraph, var_names_ntp[ipnm], var_names_ntp[opnm])
            end
        end
    end
    #* 根据有向图排序得到fluxes的计算顺序
    sorted_fluxes = AbstractComponent[]
    for idx in topological_sort(digraph)
        tmp_var_nm = var_names[idx]
        if (tmp_var_nm in output_names)
            tmp_flux = fluxes_ntp[tmp_var_nm]
            if !(tmp_flux in sorted_fluxes)
                push!(sorted_fluxes, tmp_flux)
            end
        end
    end
    sorted_fluxes
end

"""
$(SIGNATURES)

Construct a calculation graph based on all hydrological components in the hydrological unit
"""
function sort_components(components::AbstractVector{<:AbstractComponent})
    input_names, output_names, state_names = get_var_names(components)
    components_ntp = reduce(merge, map(components) do component
        tmp_input_names, tmp_output_names, tmp_state_names = get_var_names(component)
        println((tmp_input_names, tmp_output_names, tmp_state_names))
        tmp_output_state_names = vcat(tmp_output_names, tmp_state_names)
        NamedTuple{Tuple(tmp_output_state_names)}(repeat([component], length(tmp_output_state_names)))
    end)
    var_names = reduce(union, [input_names, output_names, state_names])
    var_names_ntp = namedtuple(var_names, collect(1:length(var_names)))
    digraph = SimpleDiGraph(length(var_names))
    for component in components
        tmp_input_names, tmp_output_names, tmp_state_names = get_var_names(component)
        tmp_output_names = vcat(tmp_output_names, tmp_state_names)
        for ipnm in tmp_input_names
            for opnm in tmp_output_names
                add_edge!(digraph, var_names_ntp[ipnm], var_names_ntp[opnm])
            end
        end
    end
    sorted_components = AbstractComponent[]
    for idx in topological_sort(digraph)
        tmp_var_nm = var_names[idx]
        if (tmp_var_nm in output_names)
            tmp_component = components_ntp[tmp_var_nm]
            if !(tmp_component in sorted_components)
                push!(sorted_components, tmp_component)
            end
        end
    end
    sorted_components
end

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
