"""
$(TYPEDEF)
A custom ODEProblem solver
# Fields
$(FIELDS)
"""
@kwdef struct ODESolver <: AbstractSolver
    alg::OrdinaryDiffEqAlgorithm = Tsit5()
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
    if SciMLBase.successful_retcode(sol)
        solved_states = hcat(sol.u...)
        [solved_states[i,:] for i in 1:size(solved_states)[1]]
    else
        @error "ode failed to solve"
        false
    end
end


# todo
struct DisSolver <: AbstractSolver
    
end

