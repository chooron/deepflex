function CascadeRouteFlux(
    input::Num,
)
    @parameters k [description = "水库的平均滞留时间"]
    @parameters n [description = "水库的数目"]

    return ContRouteFlux(
        input,
        [k, n],
        Num[],
        :cascade
    )
end

function (flux::ContRouteFlux{:cascade})(input::AbstractMatrix, pas::ComponentVector; kwargs...)
    n = Int(pas[:params].n)
    init_states = zeros(n)
    input_len = size(input)[2]
    input_itp = LinearInterpolation(input[1, :], collect(1:input_len))

    function nash_unithydro!(du, u, p, t)
        k = p
        du[1] = input_itp(t) - u[1] / k
        for i in 2:n
            du[i] = (u[i-1] - u[i]) / k
        end
    end

    prob = ODEProblem(nash_unithydro!, init_states, (1, length(input_vec)), (pas[:params].k,))
    sol = solve(prob, Tsit5())
    sol_vec = sol.u[:, end] .* pas[:params].k
    reshape(sol_vec, 1, input_len)
end

function get_rflux_initstates(::ContRouteFlux{:cascade}; pas::ComponentVector, ndtypes::AbstractVector{Symbol})
    reduce(vcat, [zeros(eltype(pas), Int(pas[:params][ndtype][:n])) for ndtype in ndtypes])
end

function get_rflux_func(::ContRouteFlux{:cascade}; pas::ComponentVector, ndtypes::AbstractVector{Symbol})
    #* node_num * ts_len
    node_iuh_nums = [pas[:params][ndtype][:n] for ndtype in ndtypes]
    start_idxes = Int[1; cumsum(node_iuh_nums)[1:end-1] .+ 1]
    end_idxes = Int.(cumsum(node_iuh_nums))
    iuh_states_idxes = [Int.(sidx:eidx) for (sidx, eidx) in zip(start_idxes, end_idxes)]

    function cal_q_out!(du, uh_states, q_in, q_gen, p)
        k_ps = [p[ndtype][:k] for ndtype in ndtypes]
        dstates = du[:flux_states]
        dstates[start_idxes] = @.(q_in + q_gen - uh_states[start_idxes] / k_ps)
        for (i, (k, n)) in enumerate(zip(k_ps, node_iuh_nums))
            iuh_states_idx = iuh_states_idxes[i]
            for j in 2:Int(n)
                dstates[iuh_states_idx[j]] = (uh_states[iuh_states_idx[j-1]] - uh_states[iuh_states_idx[j]]) / k
            end
        end
        q_out = [k_ps[i] * uh_states[end_idxes[i]] for i in eachindex(k_ps)]
        q_out
    end

    function cal_q_out(uh_states, q_in, q_gen, p)
        k_ps = [p[ndtype][:k] for ndtype in ndtypes]
        q_out = [k_ps[i] * uh_states[end_idxes[i]] for i in eachindex(k_ps)]
        q_out
    end

    return cal_q_out!, cal_q_out
end
