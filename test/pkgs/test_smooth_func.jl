using OrdinaryDiffEq, ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters σ ρ β
@variables x(t)

step_func(v) = (tanh(β * v) + 1.0) * 0.5

eqs = [
    D(x) ~ step_func(t - ρ) * σ
]

@mtkbuild sys = ODESystem(eqs, t)

u0 = [x => 1.0]

p = [
    σ => 1.0,
    ρ => 5.0,
    β => 1.0]

tspan = (4.5, 5.5)
prob = ODEProblem(sys, u0, tspan, p, jac=true)

for i in [0.1, 1.0, 5.0, 10.0, 50.0]
    new_prob = remake(prob, p=[β => i])
    new_sol = solve(new_prob, Tsit5(), saveat=0.1)
    println(vcat(new_sol.u...))
end