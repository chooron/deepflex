# simplified hydrological discharge model
function solve_hdm(input_vec, params)
    input_itp = LinearInterpolation(input_vec, collect(1:length(input_vec)))

    function hdm_ode!(du, u, p, t)
        s_river, q_in = u[1], u[2]
        lag = p[1]
        q_rf = (s_river + q_in) / (lag + 1)
        # no other spatial route
        q_out = q_rf + input_itp(t)
        du[1] = q_in - q_rf
        du[2] = q_out - q_in
    end

    #* init s_river and inflow
    init_states = zeros(2)
    prob = ODEProblem(hdm_ode!, init_states, (1, length(input_vec)), (params.lag,))
    sol = solve(prob, Tsit5())
    s_river_vec = sol.u[:, 1]
    q_in_vec = sol.u[:, 2]
    q_out_vec = (s_river_vec .+ q_in_vec) ./ (lag + 1) .+ input_vec
    q_out_vec
end

function get_hdm_initstates(params, ndtypes)
    [params[ndtype][:s_river] for ndtype in ndtypes]
end

function get_hdm_func(; kwargs...)
    ndtypes = kwargs[:ndtypes]
    function cal_q_out!(du, s_rivers, q_in, q_gen, p)
        lag_ps = [p[ndtype][:lag] for ndtype in ndtypes]
        q_rf = @.((s_rivers + q_in) / (lag_ps + 1))
        du[:flux_states][:] = q_in .- q_rf
        q_out = q_rf .+ q_gen
        q_out
    end

    function cal_q_out(s_rivers, q_in, q_gen, p)
        lag_ps = [p[ndtype][:lag] for ndtype in ndtypes]
        q_rf = @.((s_rivers + q_in) / (lag_ps + 1))
        q_out = q_rf .+ q_gen
        q_out
    end

    return cal_q_out!, cal_q_out
end

function DischargeRouteFlux(
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