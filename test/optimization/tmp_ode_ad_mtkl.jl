using ModelingToolkit, Optimization
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using OptimizationOptimisers
using SymbolicIndexingInterface
using SymbolicIndexingInterface: parameter_index
using SciMLStructures: replace, Tunable
using SciMLSensitivity
using Zygote

@mtkmodel FOL begin
    @parameters begin
        a
        b
    end
    @variables begin
        x(t)
        y(t)
    end
    @equations begin
        D(x) ~ (1 - x) / a + b
        D(y) ~ (1 - y) / a + b
    end
end

@mtkbuild fol = FOL()
prob = ODEProblem(fol, [fol.x => 0.0, fol.y => 0.0], (1.0, 10.0), [fol.a => 3.0, fol.b => 4.0])
sol = OrdinaryDiffEq.solve(prob, Tsit5(), saveat=1.0)
sol_u = sol.u

sol_x = hcat(sol.u...)
# ts = 1:10
# target = @.((ts - (ts^2) / 2) / 3.0 + 4.0 * ts)


# p_setter! = setp(prob, [fol.a, fol.b])
# function loss(x, p)
#     prob_ps = replace(Tunable(), SymbolicIndexingInterface.parameter_index(prob, :a), x)
#     println(prob_ps)
#     new_prob = remake(prob; p=prob_ps)
#     sol = OrdinaryDiffEq.solve(new_prob, Tsit5(), saveat=1.0)
#     sol_u = sol[1, :]
#     return sum(abs.(target .- sol_u))
# end

# optfun = OptimizationFunction(loss, Optimization.AutoForwardDiff())
# optprob = OptimizationProblem(optfun, [3.0, 4.0], ())
# solve(optprob, OptimizationOptimisers.Adam(), maxiters=10)
# SymbolicIndexingInterface.variable_index(prob, :x)
# parameter_values(prob)
# prob.ps
# replace(Tunable(), state_values(prob),[1.0])