using ModelingToolkit
using OrdinaryDiffEq
using DataInterpolations
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters α β γ δ
@variables x(t) y(t)

itp_func(t) = LinearInterpolation(rand(100), 1:100, extrapolate=true)(t)
eqs = [
    D(x) ~ (α - β * y) * x
    D(y) ~ (δ * x - γ) * y
]
odesys = structural_simplify(ODESystem(eqs, t, name=:sys))

add_eqs =
    [x ~ itp_func(t),
        y ~ itp_func(t)
    ]
compose(ODESystem(add_eqs, t; name=Symbol(ele.name, :comp_sys)), odesys)
# prob = ODEProblem(odesys, [], (0.0, 10.0), [])