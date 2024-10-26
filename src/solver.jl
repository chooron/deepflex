"""
A custom ODEProblem solver
"""
@kwdef struct ODESolver <: AbstractSolver
    alg = Tsit5()
    sensealg = InterpolatingAdjoint()
    reltol = 1e-3
    abstol = 1e-3
    saveat = 1.0
end

function (solver::ODESolver)(
    ode_func!::Function,
    pas::ComponentVector,
    initstates::AbstractArray,
    timeidx::AbstractVector;
    convert_to_array::Bool=true
)
    #* build problem
    # 虽然pas本身就包含了initstates但是initstates的构建方式因输入会有所不同
    prob = ODEProblem(
        ode_func!,
        initstates,
        (timeidx[1], timeidx[end]),
        pas
    )
    #* solve problem
    sol = solve(
        prob,
        solver.alg,
        saveat=timeidx,
        reltol=solver.reltol,
        abstol=solver.abstol,
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
@kwdef struct DiscreteSolver <: AbstractSolver
    alg = FunctionMap()
    sensealg = InterpolatingAdjoint()
end

function (solver::DiscreteSolver)(
    ode_func!::Function,
    params::ComponentVector,
    initstates::AbstractArray,
    timeidx::AbstractVector;
    convert_to_array::Bool=true
)
    #* build problem
    # 虽然pas本身就包含了initstates但是initstates的构建方式因输入会有所不同
    prob = DiscreteProblem(
        ode_func!,
        initstates,
        (timeidx[1], timeidx[end]),
        params
    )
    #* solve problem
    sol = solve(
        prob,
        solver.alg,
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

@kwdef struct ManualSolver <: AbstractSolver
    #* 计算效率过差不予考虑
end

function (solver::ManualSolver)(
    ode_func!::Function,
    pas::ComponentVector,
    initstates::AbstractVector,
    timeidx::AbstractVector
)
    T = promote_type(eltype(pas), eltype(initstates))
    init_du = zeros(T, size(initstates))
    itegration(st, pas, t) = begin
        state, states_results = st
        ode_func!(init_du, state, pas, t)
        state = state .+ init_du
        return state, (states_results..., state)
    end
    final_states, states_results = reduce((acc, t) -> itegration(acc, pas, t), timeidx, init=(initstates, (initstates,)))
    reduce(hcat, states_results)
end