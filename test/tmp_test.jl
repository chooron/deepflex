using DifferentialEquations
using ComponentArrays

function f1(u, p, t)
    @info p
    ComponentVector(a=1.01 * u[:a])
end

tspan = (0.0, 1.0)
prob = ODEProblem(f1, ComponentVector(a=1/2), tspan, (a=maximum,))
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)