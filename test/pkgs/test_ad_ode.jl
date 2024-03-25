using OrdinaryDiffEq,
    Optimization, OptimizationOptimisers, SciMLSensitivity,
    Zygote,NamedTupleTools,ComponentArrays

function lotka_volterra!(du, u, p, t)
    x, y = u
    α, β, δ, γ = p.α,p.β,p.δ,p.γ
    du[1] = dx = α * x - β * x * y
    du[2] = dy = -δ * y + γ * x * y
end

# Initial condition
u0 = [1.0, 1.0]

# Simulation interval and intermediary points
tspan = (0.0, 10.0)
tsteps = 0.0:0.1:10.0

# LV equation parameter. p = [α, β, δ, γ]
p = namedtuple([:α, :β, :δ, :γ], [1.5, 1.0, 3.0, 1.0])

# Setup the ODE problem, then solve
prob = ODEProblem(lotka_volterra!, u0, tspan, p)
sol = solve(prob, Tsit5())

# Plot the solution

function loss(p)
    sol = solve(prob, Tsit5(), p=p, saveat=tsteps)
    loss = sum(abs2, sol .- 1)
    return loss
end

callback = function (p, l)
    println(typeof(p))
    # Tell Optimization.solve to not halt the optimization. If return true, then
    # optimization stops.
    return false
end

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, ComponentVector(p))

result_ode = Optimization.solve(optprob, Adam(),
    callback=callback,
    maxiters=100)