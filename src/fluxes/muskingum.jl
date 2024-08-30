function MuskingumRouteFlux(
    input::Num,
    output::Union{Num,Nothing}=nothing,
)
    @parameters k, x, dt

    if isnothing(output)
        input_name = Symbolics.tosymbol(input, escape=false)
        output_name = Symbol(input_name, :_routed)
        output = first(@variables $output_name)
    end

    return VectorRouteFlux(
        input,
        [k, x, dt],
        routetype=:muskingum,
        output=output
    )
end

function (flux::VectorRouteFlux{:muskingum})(input::AbstractMatrix, pas::ComponentVector; kwargs...)
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

function get_rflux_initstates(::VectorRouteFlux{:muskingum}; pas::ComponentVector, ptypes::AbstractVector{Symbol})
    [pas[:initstates][ptype][:s_river] for ptype in ptypes]
end

function get_rflux_func(::VectorRouteFlux{:muskingum}; pas::ComponentVector, ptypes::AbstractVector{Symbol})
    function cal_q_out!(du, q_out, q_in, q_gen, p)
        k_ps = [p[ptype][:k] for ptype in ptypes]
        x_ps = [p[ptype][:x] for ptype in ptypes]
        dt_ps = [p[ptype][:dt] for ptype in ptypes]

        c0_ps = @.(((dt_ps / k_ps) - (2 * x_ps)) / ((2 * (1 - x_ps)) + (dt_ps / k_ps)))
        c1_ps = @.(((dt_ps / k_ps) + (2 * x_ps)) / ((2 * (1 - x_ps)) + (dt_ps / k_ps)))
        c2_ps = @.(((2 * (1 - x_ps)) - (dt_ps / k_ps)) / ((2 * (1 - x_ps)) + (dt_ps / k_ps)))

        new_q_out = @.((c0_ps * q_in) + (c1_ps * q_in[Int(t)-1]) + (c2_ps * (q_out + q_gen)))
        du[:flux_states] = new_q_out
        new_q_out
    end

    function cal_q_out(q_out, q_in, q_gen, p)
        k_ps = [p[ptype][:k] for ptype in ptypes]
        x_ps = [p[ptype][:x] for ptype in ptypes]
        dt_ps = [p[ptype][:dt] for ptype in ptypes]

        c0_ps = @.(((dt_ps / k_ps) - (2 * x_ps)) / ((2 * (1 - x_ps)) + (dt_ps / k_ps)))
        c1_ps = @.(((dt_ps / k_ps) + (2 * x_ps)) / ((2 * (1 - x_ps)) + (dt_ps / k_ps)))
        c2_ps = @.(((2 * (1 - x_ps)) - (dt_ps / k_ps)) / ((2 * (1 - x_ps)) + (dt_ps / k_ps)))

        new_q_out = @.((c0_ps * q_in) + (c1_ps * q_in[Int(t)-1]) + (c2_ps * (q_out + q_gen)))
        new_q_out
    end

    return cal_q_out!, cal_q_out
end
