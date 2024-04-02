@kwdef struct ODESolver <: AbstractSolver
    alg::OrdinaryDiffEqAlgorithm = BS3()
    # sensealg = ForwardDiffSensitivity()
    config::Dict = Dict(:reltol => 1e-3, :abstol => 1e-3)
end

function (solver::ODESolver)(
    ode_prob::ODEProblem,
    init_states::NamedTuple,
    solve_config::NamedTuple
)
    # build problem
    # tspan = (solve_config.saveat[1], solve_config.saveat[end])
    # ode_prob = ODEProblem(ode_func, init_states, tspan, func_parameters)
    # solve problem
    sol = solve(ode_prob, solver.alg,
        saveat=solve_config.saveat,
        reltol=solver.config[:reltol],
        abstol=solver.config[:abstol],
        sensealg=solver.sensealg)
    solved_u = hcat(sol.u...)
    state_names = collect(keys(init_states))
    namedtuple(state_names, [solved_u[idx, :] for idx in 1:length(state_names)])
end

@kwdef struct DiscreteSolver <: AbstractSolver
    alg = SimpleFunctionMap()
    sensealg = ForwardDiffSensitivity()
    config::Dict = Dict(:reltol => 1e-3, :abstol => 1e-3)
end

function (solver::DiscreteSolver)(
    ode_func::Function,
    func_parameters::ComponentVector{T},
    init_states::ComponentVector{T},
    solve_config::NamedTuple
) where {T<:Number}
    # build problem
    tspan = (solve_config.saveat[1], solve_config.saveat[end])
    prob = DiscreteProblem(ode_func, init_states, tspan, func_parameters)
    # solve problem
    sol = solve(prob, solver.alg,
        dt=solve_config.dt, saveat=solve_config.saveat,
        reltol=solver.config[:reltol], abstol=solver.config[:abstol],
        sensealg=solver.sensealg)
    solved_u = hcat(sol.u...)
    state_names = collect(keys(init_states))
    ComponentVector(namedtuple(state_names, [solved_u[idx, :] for idx in 1:length(state_names)]))
end

const DEFAULT_SOLVER = ODESolver()
