function solve_mskfunc(input_vec, params)
    k, x, dt = params.k, params.x, params.dt
    c0 = ((dt / k) - (2 * x)) / ((2 * (1 - x)) + (dt / k))
    c1 = ((dt / k) + (2 * x)) / ((2 * (1 - x)) + (dt / k))
    c2 = ((2 * (1 - x)) - (dt / k)) / ((2 * (1 - x)) + (dt / k))
    function msk_prob(u, p, t)
        println(t)
        q0 = u[1]
        c0, c1, c2 = p
        input1 = input_vec[Int(t)]
        input0 = input_vec[Int(t)-1]
        new_q = (c0 * input1) + (c1 * input0) + (c2 * q0)
        [new_q]
    end

    prob = DiscreteProblem(msk_prob, [input_vec[1]], (2, length(input_vec)), ComponentVector(c0=c0, c1=c1, c2=c2))
    sol = solve(prob, FunctionMap())
    reduce(vcat, sol.u)
end

function MuskingumRouteFlux(
    input::Num,
    params::Vector{Num},
)
    return RouteFlux(
        [input],
        params,
        solve_mskfunc,
        :muskingum_cunge
    )
end