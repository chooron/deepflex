struct GridRoute <: AbstractGridRoute
    "流向矩阵"
    flwdir::AbstractMatrix
    "节点位置信息"
    positions::AbstractVector
    "汇流函数"
    rfuncs::AbstractVector{<:AbstractRouteFlux}
    "bucket information: keys contains: input, output, param, state"
    infos::NamedTuple

    function GridRoute(
        name::Symbol;
        rfuncs::AbstractVector{<:AbstractRouteFlux},
        flwdir::AbtractMatrix,
        positions::AbstractVector,
    )
        #* Extract all variable names of funcs and dfuncs
        input_names, output_names, state_names = get_var_names(rfuncs)
        #* Extract all parameters names of funcs and dfuncs
        param_names = get_param_names(rfuncs)
        #* Extract all neuralnetwork names of the funcs
        nn_names = get_nn_names(rfuncs)
        #* Setup the name information of the hydrobucket
        infos = (name=name, input=input_names, output=output_names, state=state_names, param=param_names, nn=nn_names)

        return new(
            flwdir,
            positions,
            rfuncs,
            infos,
        )
    end
end

"""
input dims: node_num * ts_len
"""
function grid_routing(input::AbstractMatrix, positions::AbstractVector, flwdir::AbstractMatrix)
    #* 转换为input的稀疏矩阵
    input_arr = Array(sparse(
        [pos[1] for pos in positions],
        [pos[2] for pos in positions],
        eachslice(input, dims=1),
        size(flwdir)[1], size(flwdir)[2]
    ))
    #* 计算权重求和结果
    input_routed = sum(collect([pad_zeros(input_arr .* (flwdir .== code), arg)
                                for (code, arg) in zip(d8_codes, d8_nn_pads)]))
    #* 裁剪输入矩阵边框
    clip_arr = input_routed[2:size(input_arr)[1]+1, 2:size(input_arr)[2]+1]

    #* 将输入矩阵转换为向量
    collect([clip_arr[pos[1], pos[2]] for pos in positions])
end

function (route::GridRoute)(
    input::AbstractArray,
    pas::ComponentVector,
    timeidx::AbstractVector,
    ndtypes::AbstractVector{Symbol},
)
    @assert size(input)[1] == length(route.rfuncs) "输入变量维度需要与该模块的汇流计算维度一致"
    #* var num * node num * ts len
    #* 计算出每个node结果的插值函数
    itp_funcs = map(collect(1:size(input)[1])) do i
        LinearInterpolation.(eachslice(input[i, :, :], dims=1), Ref(timeidx), extrapolate=true)
    end
    # todo 先假设就只有一个输入
    itp_func = itp_funcs[1]

    cal_flux_q_out!, cal_flux_q_out = rfunc.qout_func(; ndtypes)
    flux_initstates = rfunc.get_initstates(; pas[:initstates], ndtypes)

    function grid_route!(du, u, p, t)
        q_out = cal_flux_q_out!(du, u[:flux_states], u[:q_in], itp_func(t), p)
        new_q_in = grid_routing(q_out, route.positions, route.flwdir)
        du[:q_in][:] = new_q_in .- q_in
    end

    init_states = ComponentVector(flux_states=flux_initstates, q_in=zeros(size(input_mat)[1]))
    prob = ODEProblem(grid_route!, init_states, (1, size(input_mat)[1]), pas[:params])
    sol = solve(prob, Tsit5())

    flux_states, q_inflows = sol[:, 1, :], sol[:, 2, :]
    q_out = cal_flux_q_out.(eachslice(flux_states, dims=2), eachslice(q_inflows, dims=2), eachslice(input, dims=2), Ref(p), Ref(true))
end

struct VectorRoute <: AbstractVectorRoute
    "流向矩阵"
    flwdir::AbstractMatrix
    "节点位置信息"
    positions::AbstractVector
    "汇流函数"
    rfuncs::AbstractVector{<:AbstractRouteFlux}

    function VectorRoute(
        name::Symbol;
        rfuncs::AbstractVector{<:AbstractRouteFlux},
        network::DiGraph,
    )
        #* Extract all variable names of funcs and dfuncs
        input_names, output_names, state_names = get_var_names(rfuncs)
        #* Extract all parameters names of funcs and dfuncs
        param_names = get_param_names(rfuncs)
        #* Extract all neuralnetwork names of the funcs
        nn_names = get_nn_names(rfuncs)
        #* Setup the name information of the hydrobucket
        infos = (name=name, input=input_names, output=output_names, state=state_names, param=param_names, nn=nn_names)

        return new(
            rfuncs,
            network,
            infos,
        )
    end
end


function (route::VectorRoute)(
    input::AbstractArray,
    pas::ComponentVector,
    timeidx::AbstractVector,
    ndtypes::AbstractVector{Symbol},
)
    #* 获取每个节点的计算顺序
    sorted_idx = topological_sort(route.network)
    input_names = [Symbol(:num_, idx) for idx in 1:size(input)[2]]
    input_ntp = NamedTuple(Tuple(input_names))(eachslice(input, dims=2))
    for cur_idx in sorted_idx
        up_idxes = inneighbors(network.network, cur_idx)
        if length(up_idxes) > 0
            #* 当前单元产流加上游单元的汇入
            cur_input = input_ntp[vcat(Symbol(:num_, cur_idx), [Symbol(:num_, up_idx) for up_idx in up_idxes])]
            #* 计算当前网格的汇出流量
            input_ntp = merge(input_ntp, route.rfuncs[1](sum(collect(cur_input)), pas[ndtypes[cur_idx]], timeidx))
        else
            #* 无上游汇入直接使用当前产流计算
            input_ntp = merge(input_ntp, route.rfuncs[1](input_ntp[Symbol(:num_, cur_idx)], pas[ndtypes[cur_idx]], timeidx))
        end
    end
    output_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), input_ntp)
    permutedims(output_arr, (1, 3, 2))
end