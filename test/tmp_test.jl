using DifferentialEquations
f(u, p, t) = 1.01 * u

u0 = (1 / 2,)
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

function f3(a::NamedTuple{T}) where {T<:Number}
    println("ok")
end