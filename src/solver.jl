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
    ode_func!::Function,
    pas::ComponentVector,
    initstates::AbstractArray,
    timeidx::AbstractVector
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
        saveat=solver.saveat,
        reltol=solver.reltol,
        abstol=solver.abstol,
        sensealg=solver.sensealg
    )
    sol_arr = Array(sol)
    if SciMLBase.successful_retcode(sol)
        return sol_arr
    else
        return false
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
    sensealg = InterpolatingAdjoint()
end

function (solver::DiscreteSolver)(
    ode_func!::Function,
    pas::ComponentVector,
    initstates::AbstractArray,
    timeidx::AbstractVector
)
    #* build problem
    # 虽然pas本身就包含了initstates但是initstates的构建方式因输入会有所不同
    prob = DiscreteProblem(
        ode_func!,
        initstates,
        (timeidx[1], timeidx[end]),
        pas
    )
    #* solve problem
    sol = solve(
        prob,
        solver.alg,
    )
    sol_arr = Array(sol)
    if SciMLBase.successful_retcode(sol)
        return sol_arr
    else
        return false
    end
end

@kwdef struct ManualSolver <: AbstractSolver
    # todo impletement the manual solver for state at each step
end

function (solver::ManualSolver)(
    ode_func!::Function,
    pas::ComponentVector,
    initstates::AbstractArray,
    timeidx::AbstractVector
)
    #* build problem
    # 虽然pas本身就包含了initstates但是initstates的构建方式因输入会有所不同
    prob = DiscreteProblem(
        ode_func!,
        initstates,
        (timeidx[1], timeidx[end]),
        pas
    )
    #* solve problem
    sol = solve(
        prob,
        solver.alg,
    )
    sol_arr = Array(sol)
    if SciMLBase.successful_retcode(sol)
        return sol_arr
    else
        return false
    end
end