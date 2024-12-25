"""
    ManualSolver{mutable} <: AbstractHydroSolver

A custom manual solver for solving ODE problems.

The `mutable` type parameter is used to indicate whether the solver uses mutable arrays (true) or immutable arrays (false).

The manual solver is a simple and lightweight solver that uses a loop to iterate over the time steps and update the state variables. It is suitable for small to medium-sized problems.

Note that setting `mutable=true` can result in a 30% performance improvement compared to `mutable=false`, since it avoids the overhead of creating new arrays at each time step.

However, it also means that the solver will modify the input arrays in-place, which may not be desirable in some cases.
"""
struct ManualSolver{mutable} end

function (solver::ManualSolver{true})(
    du_func::Function,
    pas::AbstractVector,
    initstates::AbstractArray{<:Number, 1},
    timeidx::AbstractVector;
    convert_to_array::Bool=true
)
    T1 = promote_type(eltype(pas), eltype(initstates))
    states_results = zeros(eltype(initstates), length(initstates), length(timeidx))
    tmp_initstates = copy(initstates)
    for (i, t) in enumerate(timeidx)
        tmp_du = du_func(tmp_initstates, pas, t)
        tmp_initstates = tmp_initstates .+ tmp_du
        states_results[:, i] = tmp_initstates
    end
    states_results
end

function (solver::ManualSolver{true})(
    du_func::Function,
    pas::AbstractVector,
    initstates::AbstractArray{<:Number, 2},
    timeidx::AbstractVector;
    convert_to_array::Bool=true
)
    T1 = promote_type(eltype(pas), eltype(initstates))
    states_results = zeros(eltype(initstates), size(initstates)..., length(timeidx))
    tmp_initstates = copy(initstates)
    for (i, t) in enumerate(timeidx)
        tmp_du = du_func(tmp_initstates, pas, t)
        tmp_du_mat = reduce(hcat, tmp_du)
        tmp_initstates = tmp_initstates .+ tmp_du_mat
        states_results[:, :, i] .= tmp_initstates
    end
    states_results
end

function (solver::ManualSolver{false})(
    du_func::Function,
    pas::AbstractVector,
    initstates::AbstractArray{<:Number, 1},
    timeidx::AbstractVector;
    convert_to_array::Bool=true
)
    states_results = []
    tmp_initstates = copy(initstates)
    for t in timeidx
        tmp_du = du_func(tmp_initstates, pas, t)
        tmp_initstates = tmp_initstates .+ tmp_du
        states_results = vcat(states_results, [tmp_initstates])
    end
    reduce((m1, m2) -> cat(m1, m2, dims=length(size(initstates))+1), states_results)
end

function (solver::ManualSolver{false})(
    du_func::Function,
    pas::AbstractVector,
    initstates::AbstractArray{<:Number, 2},
    timeidx::AbstractVector;
    convert_to_array::Bool=true
)
    states_results = []
    tmp_initstates = copy(initstates)
    for t in timeidx
        tmp_du = du_func(tmp_initstates, pas, t)
        tmp_initstates = tmp_initstates .+ tmp_du
        states_results = vcat(states_results, [tmp_initstates])
    end
    output = reduce((m1, m2) -> cat(m1, m2, dims=length(size(initstates))+1), states_results)
    output
end
