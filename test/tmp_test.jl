using ModelingToolkit

@variables t
D = Differential(t)

@mtkmodel FOL begin
    @parameters begin
        τ # parameters
    end
    @variables begin
        x(t) # dependent variables
    end
    @equations begin
        D(x) ~ (1 - x) / τ
    end
end

using DifferentialEquations: solve
@mtkbuild fol = FOL()
prob = ODEProblem(fol, [fol.x => 0.0], (0.0, 10.0), [fol.τ => 3.0])
sol = solve(prob)

using Plots
plot(sol)

FOL.structure[:parameters]
