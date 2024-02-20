using DifferentialEquations
using OrdinaryDiffEq
f(u, p, t) = 1.01 * u
u0 = 1 / 2
tspan = (0., 100)
prob = DiscreteProblem(f, u0, tspan)
sol = solve(prob, FunctionMap())