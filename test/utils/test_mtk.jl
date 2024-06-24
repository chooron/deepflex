using ModelingToolkit
using ModelingToolkit: t_nounits as t
using OrdinaryDiffEq: solve, FunctionMap

@inline function rate_to_proportion(r, t)
    1 - exp(-r * t)
end
@parameters c δt β γ
@constants h = 1
@variables S(t) I(t) R(t)
k = ShiftIndex(t)
infection = rate_to_proportion(
    β * c * I(k - 1) / (S(k - 1) * h + I(k - 1) + R(k - 1)), δt * h) * S(k - 1)
recovery = rate_to_proportion(γ * h, δt) * I(k - 1)

# Equations
eqs = [S(k) ~ S(k - 1) - infection * h,
    I(k) ~ I(k - 1) + infection - recovery,
    R(k) ~ R(k - 1) + recovery]
@mtkbuild sys = DiscreteSystem(eqs, t)

u0 = [S(k - 1) => 990.0, I(k - 1) => 10.0, R(k - 1) => 0.0]
p = [β => 0.05, c => 10.0, γ => 0.25, δt => 0.1]
tspan = (0.0, 100.0)
prob = DiscreteProblem(sys, u0, tspan, p)
sol = solve(prob, FunctionMap())