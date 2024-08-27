struct GridRoute <: AbstractGridRoute
    "汇流函数"
    rfunc::AbstractRouteFlux
    "流向矩阵"
    flwdir::AbstractMatrix
    "节点位置信息"
    positions::AbstractVector
    "bucket information: keys contains: input, output, param, state"
    infos::NamedTuple

    function GridRoute(
        name::Symbol;
        rfunc::AbstractRouteFlux,
        flwdir::AbstractMatrix,
        positions::AbstractVector,
    )
        #* Extract all variable names of funcs and dfuncs
        input_names, output_names, state_names = get_var_names(rfunc)
        #* Extract all parameters names of funcs and dfuncs
        param_names = get_param_names(rfunc)
        #* Extract all neuralnetwork names of the funcs
        nn_names = get_nn_names(rfunc)
        #* Setup the name information of the hydrobucket
        infos = (name=name, input=input_names, output=output_names, state=state_names, param=param_names, nn=nn_names)

        return new(
            rfunc,
            flwdir,
            positions,
            infos,
        )
    end
end

"""
input dims: node_num * ts_len
"""
function grid_routing(input::AbstractVector, positions::AbstractVector, flwdir::AbstractMatrix)
    #* 转换为input的稀疏矩阵
    input_arr = Array(sparse(
        [pos[1] for pos in positions],
        [pos[2] for pos in positions],
        input,
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
    #* var num * node num * ts len
    #* 计算出每个node结果的插值函数
    input_mat = input[1, :, :]
    itp_funcs = LinearInterpolation.(eachslice(input_mat, dims=1), Ref(timeidx), extrapolate=true)

    cal_flux_q_out!, cal_flux_q_out = get_rflux_func(route.rfunc; pas, ndtypes)
    flux_initstates = get_rflux_initstates(route.rfunc; pas, ndtypes)

    function grid_route_ode!(du, u, p, t)
        q_in = u[:q_in]
        q_gen = [itp_func(t) for itp_func in itp_funcs]
        q_out = cal_flux_q_out!(du, u[:flux_states], u[:q_in], q_gen, p)
        new_q_in = grid_routing(q_out, route.positions, route.flwdir)
        du[:q_in] = new_q_in .- q_in
    end

    init_states = ComponentVector(flux_states=flux_initstates, q_in=zeros(size(input_mat)[1]))
    prob = ODEProblem(grid_route_ode!, init_states, (1, size(input_mat)[2]), pas[:params])
    sol = solve(prob, Tsit5(), saveat=timeidx)

    # Extract flux_states and q_in for each time step
    flux_states_matrix = reduce(hcat, [u.flux_states for u in sol.u])
    q_in_matrix = reduce(hcat, [u.q_in for u in sol.u])

    q_out = cal_flux_q_out.(eachslice(flux_states_matrix, dims=2), eachslice(q_in_matrix, dims=2), eachslice(input_mat, dims=2), Ref(pas[:params]))
    q_out_mat = reduce(hcat, q_out)
    # Convert q_out_mat to 1 x mat size
    q_out_reshaped = reshape(q_out_mat, 1, size(q_out_mat)...)
    return q_out_reshaped
end

struct VectorRoute <: AbstractVectorRoute
    "流向矩阵"
    flwdir::AbstractMatrix
    "节点位置信息"
    positions::AbstractVector
    "汇流函数"
    rfunc::AbstractVector{<:AbstractRouteFlux}

    function VectorRoute(
        name::Symbol;
        rfunc::AbstractRouteFlux,
        network::DiGraph,
    )
        #* Extract all variable names of funcs and dfuncs
        input_names, output_names, state_names = get_var_names(rfunc)
        #* Extract all parameters names of funcs and dfuncs
        param_names = get_param_names(rfunc)
        #* Extract all neuralnetwork names of the funcs
        nn_names = get_nn_names(rfunc)
        #* Setup the name information of the hydrobucket
        infos = (name=name, input=input_names, output=output_names, state=state_names, param=param_names, nn=nn_names)

        return new(
            rfunc,
            network,
            infos,
        )
    end
end

"""
step route
"""
function (route::VectorRoute)(
    input::AbstractArray,
    pas::ComponentVector,
    timeidx::AbstractVector,
    ndtypes::AbstractVector{Symbol},
)
    #* var num * node num * ts len
    #* 计算出每个node结果的插值函数
    itp_func = LinearInterpolation.(eachslice(input[1, :, :], dims=1), Ref(timeidx), extrapolate=true)

    cal_flux_q_out!, cal_flux_q_out = get_rflux_func(route.rfunc; pas, ndtypes)
    flux_initstates = get_rflux_initstates(route.rfunc; pas, ndtypes)

    sorted_idxes = topological_sort(route.network)
    up_idxes = [inneighbors(network.network, cur_idx) for cur_idx in sorted_idx]

    function vec_route_ode!(du, u, p, t)
        q_in = u[:q_in]
        q_out = cal_flux_q_out!(du, u[:flux_states], u[:q_in], itp_func(t), p)
        for (sorted_idx, up_idx) in zip(sorted_idxes, up_idxes)
            q_out[sorted_idx] .+= sum(q_out[up_idx])
        end
        du[:q_in][:] = q_out .- q_in
    end

    init_states = ComponentVector(flux_states=flux_initstates, q_in=zeros(size(input_mat)[1]))
    prob = ODEProblem(grid_route_ode!, init_states, (1, size(input_mat)[1]), pas[:params])
    sol = solve(prob, Tsit5())

    flux_states, q_inflows = sol[:, 1, :], sol[:, 2, :]
    q_out = cal_flux_q_out.(eachslice(flux_states, dims=2), eachslice(q_inflows, dims=2), eachslice(input, dims=2), Ref(p), Ref(true))
    q_out
end

#= 
* entire route
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
            input_ntp = merge(input_ntp, route.rfunc[1](sum(collect(cur_input)), pas[ndtypes[cur_idx]], timeidx))
        else
            #* 无上游汇入直接使用当前产流计算
            input_ntp = merge(input_ntp, route.rfunc[1](input_ntp[Symbol(:num_, cur_idx)], pas[ndtypes[cur_idx]], timeidx))
        end
    end
    output_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), input_ntp)
    permutedims(output_arr, (1, 3, 2))
end
 =#