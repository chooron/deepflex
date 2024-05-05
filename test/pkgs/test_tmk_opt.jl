using ModelingToolkit, Optimization, ComponentArrays
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using OptimizationOptimisers
using SciMLStructures
using SciMLSensitivity

function setup_prob(
    sys::ODESystem,
    prob::ODEProblem;
    params::ComponentVector,
    init_states::ComponentVector
)
    #* setup init states
    u0 = [getproperty(sys, nm) => init_states[nm] for nm in keys(init_states)]
    u0var = ModelingToolkit.varmap_to_vars(u0, unknowns(sys))
    #* setup parameters
    p = [getproperty(sys, Symbol(nm)) => params[Symbol(nm)] for nm in ModelingToolkit.parameters(base_sys)]
    # pvar = ModelingToolkit.MTKParameters(sys, p, u0var)
    new_prob = remake(prob, p=p, u0=u0var)
    new_prob
end

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

@variables begin
    vs[1:2] = [3.0, 4.0] # , [bounds = (-2.0, 2.0)]
end
x_axes = getaxes(ComponentVector(a=1, b=2))
@parameters st=0.0

v1 = 1:10
target =  @.((v1 - (v1^2)/2)/3.0+4.0*v1)

loss = begin
    params = ComponentVector(vs, x_axes)
    
    # unew = [fol.x => st]
    # pnew = [fol.a => params[:a], fol.b => params[:b]]
    # unew = ModelingToolkit.varmap_to_vars([fol.x => st], unknowns(fol))
    # pnew = ModelingToolkit.MTKParameters(sys, [fol.a => params[:a], fol.b => params[:b]], unew)
    new_prob = remake(prob, u0=[fol.x => st], p=[fol.a => params[:a], fol.b => params[:b]])
    sol = OrdinaryDiffEq.solve(new_prob, Tsit5(), saveat=1.0, sensealg =ForwardDiffSensitivity())
    sol_u = [v[1] for v in sol.u]
    sum(abs.(target .- sol_u))
end

@mtkbuild sys = OptimizationSystem(loss, vs, [st]) # , constraints=cons
u0 = [vs[1] => 3.0
    vs[2] => 4.0]
prob = OptimizationProblem(sys,
    u0,
    grad=true,
    hess=true,
    cons_j=true,
    cons_h=true)
solve(prob, OptimizationOptimisers.Adam(), maxiters=10)

