
default_callback_func(p, l) = begin
    @info l
    false
end

#* build predict function
function predict_func(x::AbstractVector{T}, p) where {T}
    component, input, timeidx, solver, tunable_pas_axes, default_model_pas = p
    tmp_tunable_pas = ComponentVector(x, tunable_pas_axes)
    tmp_pas = merge_ca(default_model_pas, tmp_tunable_pas)
    component(input, tmp_pas, timeidx=timeidx, solver=solver)
end


"""
$(SIGNATURES)

Parameter optimization for global search of hydrological units, nodes
"""
function param_box_optim(
    component::AbstractComponent;
    tunable_pas::ComponentArray,
    const_pas::ComponentArray,
    input::NamedTuple,
    target::NamedTuple,
    timeidx::Vector,
    kwargs...,
)
    """
    针对模型参数进行二次优化
    """
    # 获取需要优化的参数名称
    solve_alg = get(kwargs, :solve_alg, BBO_adaptive_de_rand_1_bin_radiuslimited())
    target_name = get(kwargs, :target_name, :flow)
    loss_func = get(kwargs, :loss_func, mse)
    callback_func = get(kwargs, :callback_func, default_callback_func)
    lb = get(kwargs, :lb, zeros(length(tunable_pas)))
    ub = get(kwargs, :ub, ones(length(tunable_pas)) .* 100)
    maxiters = get(kwargs, :maxiters, 10)
    solver = get(kwargs, :solver, ODESolver(saveat=timeidx))

    default_model_pas = ComponentArray(merge_recursive(NamedTuple(tunable_pas), NamedTuple(const_pas)))
    tunable_pas_axes = getaxes(tunable_pas)

    function objective(x::AbstractVector{T}, p) where {T}
        loss = loss_func(target[target_name], getproperty(predict_func(x, p), target_name))
        loss
    end

    # 构建问题
    optf = Optimization.OptimizationFunction(objective)
    prob_args = (component, input, timeidx, solver, tunable_pas_axes, default_model_pas)
    optprob = Optimization.OptimizationProblem(optf, collect(tunable_pas), prob_args, lb=lb, ub=ub)
    sol = Optimization.solve(optprob, solve_alg, callback=callback_func, maxiters=maxiters)

    ComponentVector(sol.u, tunable_pas_axes)
end

"""
$(SIGNATURES)

Parameter optimization for local search of hydrological units, nodes
"""
function param_grad_optim(
    component::AbstractComponent;
    tunable_pas::AbstractVector,
    const_pas::AbstractVector,
    input::NamedTuple,
    target::NamedTuple,
    timeidx::Vector,
    kwargs...,
)
    """
    using Optimization.jl for param optimize
    """
    # 获取需要优化的参数名称
    solve_alg = get(kwargs, :solve_alg, Adam())
    adtype = get(kwargs, :adtype, Optimization.AutoForwardDiff())
    loss_func = get(kwargs, :loss_func, mse)
    callback_func = get(kwargs, :callback_func, default_callback_func)
    maxiters = get(kwargs, :maxiters, 10)
    solver = get(kwargs, :solver, ODESolver())

    default_model_pas = ComponentArray(merge_recursive(NamedTuple(tunable_pas), NamedTuple(const_pas)))
    tunable_pas_axes = getaxes(tunable_pas)

    function objective(x::AbstractVector{T}, p) where {T}
        predict_result = predict_func(x, p)
        loss = mean([loss_func(target[key], predict_result[key]) for key in keys(target)])
        loss
    end
    # build optim problem
    optf = Optimization.OptimizationFunction(objective, adtype)
    prob_args = (component, input, timeidx, solver, tunable_pas_axes, default_model_pas)
    optprob = Optimization.OptimizationProblem(optf, collect(tunable_pas), prob_args)
    sol = Optimization.solve(optprob, solve_alg, callback=callback_func, maxiters=maxiters)

    ComponentVector(sol.u, tunable_pas_axes)
end

function nn_param_optim(
    flux::AbstractNeuralFlux;
    input::NamedTuple,
    target::NamedTuple,
    init_params::ComponentVector,
    kwargs...
)
    solve_alg = get(kwargs, :solve_alg, Adam(0.01))
    adtype = get(kwargs, :adtype, Optimization.AutoZygote())
    maxiters = get(kwargs, :maxiters, 100)
    loss_func = get(kwargs, :loss_func, mse)
    callback_func = get(kwargs, :callback_func, default_callback_func)
    
    x = reduce(hcat, [input[k] for k in get_input_names(flux)])

    function pred_func(u)
        NamedTuple{Tuple(get_output_names(flux))}(flux(x, u))
    end

    function objective(u, p)
        predict = pred_func(u)
        mean([loss_func(predict[nm], target[nm]) for nm in keys(target)])
    end

    optf = Optimization.OptimizationFunction(objective, adtype)
    optprob = Optimization.OptimizationProblem(optf, init_params)
    sol = Optimization.solve(optprob, solve_alg,callback=callback_func, maxiters=maxiters)
    sol.u
end