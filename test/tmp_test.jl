using OrdinaryDiffEq

function f(u, p, t)
    @info t
    1.01 * u
end
u0 = 1 / 2
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8,saveat=0.0:0.001:1.0)