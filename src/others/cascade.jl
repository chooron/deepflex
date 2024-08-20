function solve_nashuh(input_vec, params)
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

    prob = ODEProblem(nash_unithydro!, init_states, (1, length(input_vec)), (params.p,))
    sol = solve(prob, Tsit5())
    sol.u
end

function CascadeRoute(
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