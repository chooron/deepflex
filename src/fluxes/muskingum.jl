function MuskingumRouteFlux(
    input::Num,
    output::Union{Num,Nothing}=nothing,
)
    @parameters k, x

    if isnothing(output)
        input_name = Symbolics.tosymbol(input, escape=false)
        output_name = Symbol(input_name, :_routed)
        output = first(@variables $output_name)
    end

    return RouteFlux(
        input,
        [k, x],
        Num[],
        routetype=:muskingum,
        output=output
    )
end

function (flux::RouteFlux{:muskingum})(input::Matrix, pas::ComponentVector, timeidx::AbstractVector; kwargs...)
    input_len = size(input)[2]
    input_itp = LinearInterpolation(input[1, :], collect(1:input_len))
    params = pas[:params]

    function msk_prob!(du, u, p, t)
        s_river = u[1]
        q_in = input_itp(t)
        k, x = p
        q_out = (s_river - k * x * q_in) / (k * (1 - x))
        du[1] = q_in - q_out
    end

    init_states = [params.k * input[1, 1]]
    prob = ODEProblem(msk_prob!, init_states, (1, input_len), params)
    sol = solve(prob, Rosenbrock23(), saveat=timeidx)

    s_river_vec = Array(sol)
    q_out_vec = @.((s_river_vec - params.k * params.x * input) / (params.k * (1 - params.x)))
    q_out_vec
end

function get_rflux_initstates(::RouteFlux{:muskingum}; input::AbstractMatrix, pas::ComponentVector, ptypes::AbstractVector{Symbol})
    [pas[:params][ptype][:k] for ptype in ptypes] .* input[:, 1]
end

function get_rflux_func(::RouteFlux{:muskingum}; pas::ComponentVector, ptypes::AbstractVector{Symbol})

    function cal_q_out!(du, s_river, q_in, q_gen, p)
        k_ps = [p[ptype][:k] for ptype in ptypes]
        x_ps = [p[ptype][:x] for ptype in ptypes]

        q_out = @.((s_river - k_ps * x_ps * q_in) / (k_ps * (1 - x_ps)))
        du[:flux_states] = q_in .- q_out
        q_out .+ q_gen
    end

    function cal_q_out(s_river, q_in, q_gen, p)
        k_ps = [p[ptype][:k] for ptype in ptypes]
        x_ps = [p[ptype][:x] for ptype in ptypes]

        q_out = @.((s_river - k_ps * x_ps * q_in) / (k_ps * (1 - x_ps)))
        q_out .+ q_gen
    end

    return cal_q_out!, cal_q_out
end
