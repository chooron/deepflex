# simplified hydrological discharge model
#= 
Hydrological discharge model: (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023WR036170)
=#
function DischargeRouteFlux(
    input::Num,
    output::Union{Num,Nothing}=nothing,
)
    @parameters lag
    @variables s_river

    if isnothing(output)
        input_name = Symbolics.tosymbol(input, escape=false)
        output_name = Symbol(input_name, :_routed)
        output = first(@variables $output_name)
    end

    return RouteFlux(
        input,
        [lag],
        [s_river],
        routetype=:discharge,
        output=output
    )
end

function (flux::RouteFlux{:discharge})(input::AbstractMatrix, pas::ComponentVector; kwargs...)
    input_len = size(input)[2]
    input_itp = LinearInterpolation(input[1, :], collect(1:input_len))

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
    reshape(q_out_vec, 1, input_len)
end

function get_rflux_initstates(::RouteFlux{:discharge}; input::AbstractMatrix, pas::ComponentVector, ptypes::AbstractVector{Symbol})
    [pas[:initstates][ptype][:s_river] for ptype in ptypes]
end

function get_rflux_func(::RouteFlux{:discharge}; pas::ComponentVector, ptypes::AbstractVector{Symbol})
    function cal_q_out!(du, s_rivers, q_in, q_gen, p)
        lag_ps = [p[ptype][:lag] for ptype in ptypes]
        q_rf = @.((s_rivers + q_in) / (lag_ps + 1))
        du[:flux_states] = q_in .- q_rf
        q_out = q_rf .+ q_gen
        q_out
    end

    function cal_q_out(s_rivers, q_in, q_gen, p)
        lag_ps = [p[ptype][:lag] for ptype in ptypes]
        q_rf = @.((s_rivers + q_in) / (lag_ps + 1))
        q_out = q_rf .+ q_gen
        q_out
    end

    return cal_q_out!, cal_q_out
end
