"""
A custom ODEProblem solver
"""
@kwdef struct ODESolver <: AbstractHydroSolver
    alg = Tsit5()
    sensealg = InterpolatingAdjoint()
    reltol = 1e-3
    abstol = 1e-3
    saveat = 1.0
end

function (solver::ODESolver)(
    du_func::Function,
    pas::ComponentVector,
    initstates::AbstractArray,
    timeidx::AbstractVector;
    convert_to_array::Bool=true
)
    ode_func! = (du, u, p, t) -> (du[:] = du_func(u, p, t))

    #* build problem
    prob = ODEProblem(ode_func!, initstates, (timeidx[1], timeidx[end]), pas)
    #* solve problem
    sol = solve(
        prob, solver.alg, saveat=timeidx,
        reltol=solver.reltol, abstol=solver.abstol,
        sensealg=solver.sensealg
    )
    if convert_to_array
        if SciMLBase.successful_retcode(sol)
            sol_arr = Array(sol)
        else
            @warn "ODE solver failed, please check the parameters and initial states, or the solver settings"
            sol_arr = zeros(size(initstates)..., length(timeidx))
        end
        return sol_arr
    else
        return sol
    end
end

"""
A custom ODEProblem solver
"""
@kwdef struct DiscreteSolver <: AbstractHydroSolver
    alg = FunctionMap{true}()
    sensealg = InterpolatingAdjoint()
end

function (solver::DiscreteSolver)(
    du_func::Function,
    params::ComponentVector,
    initstates::AbstractArray,
    timeidx::AbstractVector;
    convert_to_array::Bool=true
)
    ode_func! = (du, u, p, t) -> (du[:] = du_func(u, p, t))
    #* build problem
    prob = DiscreteProblem(ode_func!, initstates, (timeidx[1], timeidx[end]), params)
    #* solve problem
    sol = solve(prob, solver.alg, saveat=timeidx, sensealg=solver.sensealg)
    if convert_to_array
        if SciMLBase.successful_retcode(sol)
            sol_arr = Array(sol)
        else
            @warn "ODE solver failed, please check the parameters and initial states, or the solver settings"
            sol_arr = zeros(size(initstates)..., length(timeidx))
        end
        return sol_arr
    else
        return sol
    end
end

"""
    ManualSolver{mutable} <: AbstractHydroSolver

A custom manual solver for solving ODE problems.

The `mutable` type parameter is used to indicate whether the solver uses mutable arrays (true) or immutable arrays (false).

The manual solver is a simple and lightweight solver that uses a loop to iterate over the time steps and update the state variables. It is suitable for small to medium-sized problems.

Note that setting `mutable=true` can result in a 30% performance improvement compared to `mutable=false`, since it avoids the overhead of creating new arrays at each time step.

However, it also means that the solver will modify the input arrays in-place, which may not be desirable in some cases.
"""
struct ManualSolver{mutable} <: AbstractHydroSolver end

function (solver::ManualSolver{true})(
    du_func::Function,
    pas::ComponentVector,
    initstates::AbstractArray,
    timeidx::AbstractVector
)
    T1 = promote_type(eltype(pas), eltype(initstates))
    states_results = zeros(T1, size(initstates)..., length(timeidx))
    tmp_initstates = copy(initstates)
    for (i, t) in enumerate(timeidx)
        tmp_du = du_func(tmp_initstates, pas, t)
        tmp_initstates = tmp_initstates .+ tmp_du
        states_results[:, i] = tmp_initstates
    end
    states_results
end

function (solver::ManualSolver{false})(
    du_func::Function,
    pas::ComponentVector,
    initstates::AbstractArray,
    timeidx::AbstractVector
)
    states_results = []
    tmp_initstates = copy(initstates)
    for t in timeidx
        tmp_du = du_func(tmp_initstates, pas, t)
        tmp_initstates = tmp_initstates .+ tmp_du
        states_results = vcat(states_results, tmp_initstates)
    end
    reduce((m1, m2) -> cat(m1, m2, dims=length(size(initstates))+1), states_results)
end
