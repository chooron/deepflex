@kwdef struct ODESolver <: AbstractSolver
    dt::Number
    saveat::Vector
    alg::OrdinaryDiffEqAlgorithm = BS3()
    sensealg = ForwardDiffSensitivity()
    config::Dict = Dict(:reltol => 1e-3, :abstol => 1e-3)
end

function (solver::ODESolver)(ode_func::Function,
    ode_parameters::NamedTuple,
    init_states::ComponentVector{T}) where {T<:Number}
    # build problem
    tspan = (solver.saveat[1], solver.saveat[end])
    prob = ODEProblem(ode_func, init_states, tspan, ode_parameters)
    # solve problem
    sol = solve(prob, solver.alg,
        dt=solver.dt, saveat=solver.saveat,
        reltol=solver.config[:reltol], abstol=solver.config[:abstol],
        sensealg=solver.sensealg)
    solved_u = hcat(sol.u...)
    state_names = collect(keys(init_states))
    ComponentVector(namedtuple(state_names, [solved_u[idx, :] for idx in 1:length(state_names)]))
end

@kwdef struct DiscreteSolver <: AbstractSolver
    dt::Number
    saveat::Vector
    alg = SimpleFunctionMap()
    sensealg = ForwardDiffSensitivity()
    config::Dict = Dict(:reltol => 1e-3, :abstol => 1e-3)
end

function (solver::ODESolver)(ode_func::Function,
    ode_parameters::NamedTuple,
    init_states::ComponentVector{T}) where {T<:Number}
    # build problem
    tspan = (solver.saveat[1], solver.saveat[end])
    prob = DiscreteProblem(ode_func, init_states, tspan, ode_parameters)
    # solve problem
    sol = solve(prob, solver.alg,
        dt=solver.dt, saveat=solver.saveat,
        reltol=solver.config[:reltol], abstol=solver.config[:abstol])
    solved_u = hcat(sol.u...)
    state_names = collect(keys(init_states))
    ComponentVector(namedtuple(state_names, [solved_u[idx, :] for idx in 1:length(state_names)]))
end
