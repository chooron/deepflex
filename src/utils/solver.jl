@kwdef struct ODESolver <: AbstractSolver
    alg::OrdinaryDiffEqAlgorithm = Tsit5()
    sensealg = ForwardDiffSensitivity()
    reltol = 1e-3
    abstol = 1e-3
    saveat = 1.0
end

function (solver::ODESolver)(
    ode_prob::ODEProblem,
    state_names::AbstractVector{Symbol},
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
        @info "ode solved successfully"
    else
        @error "ode failed to solve"
    end
    solved_u = hcat(sol.u...)
    ComponentVector(namedtuple(state_names, [solved_u[idx, :] for idx in 1:length(state_names)]))
end


