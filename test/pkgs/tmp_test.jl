using ModelingToolkit
using OrdinaryDiffEq
using DataInterpolations
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters α β γ δ
@variables x(t) y(t)

itp_func(t) = LinearInterpolation(rand(100), 1:100, extrapolate=true)(t)
eqs = [
    D(x) ~ (α - β * y) * x + itp_func(t)
    D(y) ~ (δ * x - γ) * y + itp_func(t)
]
@mtkbuild odesys = ODESystem(eqs, t)

function f1(a, b; kwargs...)
    return a + b + get(kwargs, :c, 0)
end

f1(1,2,d=2)

struct A1
    a
end

