function solve_hdm(input_vec, params)
    n = params.n
    init_states = zeros(n)
    input_itp = LinearInterpolation(input_vec, collect(1:length(input_vec)))

    function nash_unithydro!(du, u, p, t)
        k = p
        du[1] = input_itp(t) - u[1] / k
        for i in 2:n
            du[i] = (u[i-1] - u[i]) / k
        end
    end

    prob = ODEProblem(nash_unithydro!, init_states, (1, length(input_vec)), (params.k,))
    sol = solve(prob, Tsit5())
    sol.u[:, end] .* params.k
end


function get_nashiuh_func(param_vec::AbstractVector)
    #* node_num * ts_len
    node_iuh_nums = [params[:n] for params in param_vec]
    start_idxes = [1; cumsum(node_iuh_nums)[1:end-1] .+ 1]
    end_idxes = cumsum(node_iuh_nums)
    iuh_states_idxes = [sidx:eidx for (sidx, eidx) in zip(start_idxes, end_idxes)]

    function cal_q_out!(du, uh_states, q_in, q_gen, p)
        k_ps = [p[ndtype][:k] for ndtype in ndtypes]
        dstates = du[:uh_states]
        for i in 1:length(param_vec)
            iuh_states_idx = iuh_states_idxes[i]
            dstates[iuh_states_idx[1]] = q_in[i] .+ q_gen - uh_states[iuh_states_idx[1]] / k_ps[i]
            for j in 2:node_iuh_nums[i]
                dstates[iuh_states_idx[j]] = (uh_states[iuh_states_idx[j-1]] - uh_states[iuh_states_idx[j]]) / k_ps[i]
            end
        end
        q_out = [k_ps[i] * uh_states[end_idxes[i]] for i in 1:length(param_vec)]
        q_out
    end

    function cal_q_out(uh_states, q_in, q_gen, p)
        k_ps = [p[ndtype][:k] for ndtype in ndtypes]
        q_out = [k_ps[i] * uh_states[end_idxes[i]] for i in 1:length(param_vec)]
        q_out
    end

    return cal_q_out!, cal_q_out
end

function CascadeRouteFlux(
    inputs::Num,
)
    @parameters k [description = "水库的平均滞留时间"]
    @parameters n [description = "水库的数目"]

    return RouteFlux(
        [inputs],
        [k, n],
        solve_nashuh,
        :nash_cascade
    )
end