using DifferentialEquations
using ComponentArrays

function f!(du, u, p, t)
    du = 1.01 * u[:a]
    println("du:$(du), u:$(u)")
    ComponentVector(a=du)
end

u0 = ComponentVector(a=1 / 2)
tspan = (0.0, 1.0)
prob = ODEProblem(f!, u0, tspan)
sol = solve(prob)