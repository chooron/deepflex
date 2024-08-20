
function solve_nashuh(input_vec, params, timeidx)
    n = params.n
    init_states = zeros(n)
    input_itp = LinearInterpolation(input_vec, timeidx)

    function nash_unithydro!(du, u, p, t)
        k = p
        du[1] = input_itp(t) - u[1] / k
        for i in 2:n
            du[i] = (u[i-1] - u[i]) / k
        end
    end

    prob = ODEProblem(nash_unithydro!, init_states, (timeidx[1], timeidx[end]), (params.p,))
    sol = solve(prob, Tsit5())
    sol.u
end