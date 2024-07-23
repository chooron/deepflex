"""
$(TYPEDEF)
A custom ODEProblem solver
# Fields
$(FIELDS)
"""
@kwdef struct ODESolver <: AbstractSolver
    alg = Tsit5()
    sensealg = InterpolatingAdjoint()
    reltol = 1e-3
    abstol = 1e-3
    saveat = 1.0
end

function (solver::ODESolver)(
    ode_prob::ODEProblem,
)
    sol = solve(
        ode_prob,
        solver.alg,
        saveat=solver.saveat,
        reltol=solver.reltol,
        abstol=solver.abstol,
        sensealg=solver.sensealg
    )
    num_u = length(ode_prob.u0)
    if SciMLBase.successful_retcode(sol)
        return [sol[i, :] for i in 1:num_u]
    else
        println("solve fail")
        return [zeros(length(solver.saveat)) for _ in 1:num_u]
    end
end

"""
$(TYPEDEF)
A custom ODEProblem solver
# Fields
$(FIELDS)
"""
@kwdef struct DiscreteSolver <: AbstractSolver
    alg = FunctionMap()
    reltol = 1e-3
    abstol = 1e-3
    sensealg = InterpolatingAdjoint()
end

function (solver::DiscreteSolver)(
    ode_prob::DiscreteProblem
)
    sol = solve(
        ode_prob,
        solver.alg,
        reltol=solver.reltol,
        abstol=solver.abstol,
        sensealg=solver.sensealg
    )
    num_u = length(ode_prob.u0)
    [sol[i, :] for i in 1:num_u]
end

@kwdef struct ManualSolver <: AbstractSolver
    # todo impletement the manual solver for state at each step
end