function MuskingumRouteFlux(
    input::Num,
)
    @parameters k, x, dt

    return DiscRouteFlux(
        input,
        [k, x, dt],
        Num[],
        :muskingum
    )
end

function (flux::DiscRouteFlux{:muskingum})(input::AbstractMatrix, pas::ComponentVector; kwargs...)
    input_len = size(input)[2]
    input_vec = input[1, :]
    params = pas[:params]

    k, x, dt = params.k, params.x, params.dt
    c0 = ((dt / k) - (2 * x)) / ((2 * (1 - x)) + (dt / k))
    c1 = ((dt / k) + (2 * x)) / ((2 * (1 - x)) + (dt / k))
    c2 = ((2 * (1 - x)) - (dt / k)) / ((2 * (1 - x)) + (dt / k))

    function msk_prob(u, p, t)
        q_in, q_out = u[1], u[2]
        c0, c1, c2 = p
        new_q = (c0 * input_vec[Int(t)]) + (c1 * q_in) + (c2 * q_out)
        [input_vec[Int(t)], new_q]
    end

    prob = DiscreteProblem(msk_prob, [input_vec[1], input_vec[1]], (1, input_len), ComponentVector(c0=c0, c1=c1, c2=c2))
    sol = solve(prob, FunctionMap())
    Array(sol[2,:])
end

function get_rflux_initstates(::DiscRouteFlux{:muskingum}; pas::ComponentVector, ndtypes::AbstractVector{Symbol})
    [pas[:initstates][ndtype][:s_river] for ndtype in ndtypes]
end

function get_rflux_func(::DiscRouteFlux{:muskingum}; pas::ComponentVector, ndtypes::AbstractVector{Symbol})
    function cal_q_out!(du, q_out, q_in, q_gen, p)
        k_ps = [p[ndtype][:k] for ndtype in ndtypes]
        x_ps = [p[ndtype][:x] for ndtype in ndtypes]
        dt_ps = [p[ndtype][:dt] for ndtype in ndtypes]

        c0_ps = @.(((dt_ps / k_ps) - (2 * x_ps)) / ((2 * (1 - x_ps)) + (dt_ps / k_ps)))
        c1_ps = @.(((dt_ps / k_ps) + (2 * x_ps)) / ((2 * (1 - x_ps)) + (dt_ps / k_ps)))
        c2_ps = @.(((2 * (1 - x_ps)) - (dt_ps / k_ps)) / ((2 * (1 - x_ps)) + (dt_ps / k_ps)))

        new_q_out = @.((c0_ps * q_in) + (c1_ps * q_in[Int(t)-1]) + (c2_ps * (q_out + q_gen)))
        du[:flux_states] = new_q_out
        new_q_out
    end

    function cal_q_out(q_out, q_in, q_gen, p)
        k_ps = [p[ndtype][:k] for ndtype in ndtypes]
        x_ps = [p[ndtype][:x] for ndtype in ndtypes]
        dt_ps = [p[ndtype][:dt] for ndtype in ndtypes]

        c0_ps = @.(((dt_ps / k_ps) - (2 * x_ps)) / ((2 * (1 - x_ps)) + (dt_ps / k_ps)))
        c1_ps = @.(((dt_ps / k_ps) + (2 * x_ps)) / ((2 * (1 - x_ps)) + (dt_ps / k_ps)))
        c2_ps = @.(((2 * (1 - x_ps)) - (dt_ps / k_ps)) / ((2 * (1 - x_ps)) + (dt_ps / k_ps)))

        new_q_out = @.((c0_ps * q_in) + (c1_ps * q_in[Int(t)-1]) + (c2_ps * (q_out + q_gen)))
        new_q_out
    end

    return cal_q_out!, cal_q_out
end
