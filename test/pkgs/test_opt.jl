using ModelingToolkit
using OrdinaryDiffEq
# using DataInterpolations
using Interpolations
using ModelingToolkit: t_nounits as t, D_nounits as D


@parameters α β γ δ
@variables x(t) y(t)
itp_func = t -> (sin(t) + 1) * 0.5


eqs = [
    D(x) ~ (α - β * y) * x + itp_func(t)
    D(y) ~ (δ * x - γ) * y + itp_func(t)
]
@mtkbuild odesys = ODESystem(eqs, t)

odeprob = ODEProblem(odesys, [x => 1.0, y => 1.0], (0.0, 10.0), [α => 1.5, β => 1.0, γ => 3.0, δ => 1.0])
timesteps = 0.0:0.1:10.0
sol = solve(odeprob, Tsit5(); saveat=timesteps)

data = Array(sol)
# add some random noise
data = data + 0.01 * randn(size(data))

using SymbolicIndexingInterface: parameter_values, state_values
using SciMLStructures: Tunable, replace, replace!

function loss(x, p)
    odeprob = p[1] # ODEProblem stored as parameters to avoid using global variables
    ps = parameter_values(odeprob) # obtain the parameter object from the problem
    ps = replace(Tunable(), ps, x) # create a copy with the values passed to the loss function
    # remake the problem, passing in our new parameter object
    newprob = remake(odeprob; p=ps)
    timesteps = p[2]
    sol = solve(newprob, AutoTsit5(Rosenbrock23()); saveat=timesteps)
    truth = p[3]
    data = Array(sol)
    return sum((truth .- data) .^ 2) / length(truth)
end

using Optimization
using OptimizationOptimJL

# manually create an OptimizationFunction to ensure usage of `ForwardDiff`, which will
# require changing the types of parameters from `Float64` to `ForwardDiff.Dual`
optfn = OptimizationFunction(loss, Optimization.AutoForwardDiff())
# parameter object is a tuple, to store differently typed objects together
optprob = OptimizationProblem(
    optfn, rand(4), (odeprob, timesteps, data), lb=0.1zeros(4), ub=3ones(4)) # 
sol = solve(optprob, BFGS())