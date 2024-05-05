using ModelingToolkit, Optimization, ComponentArrays
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using OptimizationOptimisers

@mtkmodel FOL begin
    @parameters begin
        a
        b # parameters
    end
    @variables begin
        x(t) # dependent variables
    end
    @equations begin
        D(x) ~ (1 - x) / a + b
    end
end

@mtkbuild fol = FOL()
prob = ODEProblem(fol, [fol.x => 0.0], (1.0, 10.0), [fol.a => 3.0, fol.b => 4.0])
sol = OrdinaryDiffEq.solve(prob, Tsit5(), saveat=1.0)

begin
    vs = [3.0, 4.0] # , [bounds = (-2.0, 2.0)]
end
x_axes = getaxes(ComponentVector(a=1, b=2))
st=0.0

v1 = 1:10
target =  @.((v1 - (v1^2)/2)/3.0+4.0*v1)
function loss(vs, st)
    params = ComponentVector(vs, x_axes)
    unew = ModelingToolkit.varmap_to_vars([fol.x => st], unknowns(fol))
    pnew = ModelingToolkit.MTKParameters(fol, [fol.a => params[:a], fol.b => params[:b]], unew)
    new_prob = remake(prob, u0=unew, p=pnew)
    sol = OrdinaryDiffEq.solve(new_prob, Tsit5(), saveat=1.0)
    sol_u = sol[1,:]
    return sum(abs.(target .- sol_u))
end

u0 = [3.0, 4.0]
optfun = OptimizationFunction(loss, Optimization.AutoForwardDiff())
optprob = OptimizationProblem(optfun,
    u0, 
    st)
solve(optprob, OptimizationOptimisers.Adam(), maxiters=10)