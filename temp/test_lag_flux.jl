using OrdinaryDiffEq

ifelse_func(x) = ifelse(x > 0, 1, 0)

uh_1_half(input, lag) = (ifelse_func(input - lag) + ifelse_func(lag - input) * ifelse_func(input) * (input / lag)^2.5)

uh_2_full(input, lag) = begin
    half_lag = lag ./ 2
    (ifelse_func(input - lag) * 1 +
     ifelse_func(lag - input) * ifelse_func(input - half_lag) * (1 - 0.5 * abs(2 - input / half_lag)^2.5) +
     ifelse_func(half_lag - input) * ifelse_func(input) * (0.5 * abs(input / half_lag)^2.5))
end

function flux_unit_routing(flux::Vector, timeidx::Vector, lag_func::Function, lag_time::T, num_uh::Int) where {T}
    ts = 1:num_uh
    lag_weights = [lag_func(t, lag_time) for t in ts]
    lag_weights = vcat([lag_weights[1]], (circshift(lag_weights, -1).-lag_weights)[1:end-1])

    function lag_disc_func(u, p, t)
        u = circshift(u, -1)
        u[end] = 0.0
        tmp_u = flux[Int(t)] .* p[:weights] .+ u
        tmp_u
    end

    prob = DiscreteProblem(lag_disc_func, lag_weights, (timeidx[1], timeidx[end]), ComponentVector(weights=lag_weights))
    sol = solve(prob, FunctionMap())
    if !SciMLBase.successful_retcode(sol)
        println("Error in GR4J routing ")
        return flux
    end
    # sol[1, :]
    sol
end

re = flux_unit_routing([1, 2, 3, 4, 2, 3, 2, 4, 2, 1], collect(1:10), uh_1_half, 3.5, 10)

f(u, p, t) = 1.01 * u
u0 = 1 / 2
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)

function ff(u, p, t)
    u1, u2 = u[:u1], u[:u2]
    ComponentVector(u1=u1 .* 1.02, u1=u2 ./ 1.02)
end

tspan = (0.0, 1.0)
prob = ODEProblem(f, ComponentVector(u1=[2 3 4;], u2=[2 3;]), tspan)
solve(prob, Tsit5())