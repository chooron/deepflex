mse(y, y_hat) = sum((y .- y_hat) .^ 2) / length(y)
default_callback_func(p, l) = begin
    @info l
    false
end
using ForwardDiff

function param_box_optim(
    component::AbstractComponent;
    tunable_pas::ComponentArray,
    const_pas::ComponentArray,
    input::NamedTuple,
    target::NamedTuple,
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

    tunable_pas_axes = getaxes(tunable_pas)

    predict_func(x, p) = begin
        tunable_pas_type = eltype(x)
        tmp_tunable_pas = ComponentVector(x, tunable_pas_axes)
        const_pas = tunable_pas_type.(const_pas)
        tmp_pas = merge_ca(tmp_tunable_pas, const_pas)[:param]
        component(input, tmp_pas)
    end

    objective(x, p) = loss_func(target[target_name], predict_func(x, p)[target_name])

    # 构建问题
    optf = Optimization.OptimizationFunction(objective)
    optprob = Optimization.OptimizationProblem(optf, collect(tunable_pas), lb=lb, ub=ub)
    sol = Optimization.solve(optprob, solve_alg, callback=callback_func, maxiters=maxiters)

    ComponentVector(sol.u, tunable_pas_axes)
end

function param_grad_optim(
    component::AbstractComponent;
    tunable_pas::AbstractVector,
    const_pas::AbstractVector,
    input::ComponentVector,
    target::ComponentVector,
    kwargs...,
)
    """
    using Optimization.jl for param optimize
    """
    # 获取需要优化的参数名称
    solve_alg = get(kwargs, :solve_alg, Adam())
    adtype = get(kwargs, :adtype, Optimization.AutoForwardDiff())  # AutoForwardDiff and AutoFiniteDiff AutoZygote
    target_name = get(kwargs, :target_name, :flow)
    loss_func = get(kwargs, :loss_func, mse)
    callback_func = get(kwargs, :callback_func, default_callback_func)
    maxiters = get(kwargs, :maxiters, 10)

    tunable_pas_axes = getaxes(tunable_pas)

    # build predict function
    function predict_func(x::AbstractVector{T}, p) where {T}
        tmp_tunable_pas = ComponentVector(x, tunable_pas_axes)
        tmp_pas = merge_ca(tmp_tunable_pas, T.(const_pas))[:param]
        component(T.(input), tmp_pas)
    end

    function objective(x::AbstractVector{T}, p) where {T}
        loss_func(T.(target[target_name]), predict_func(x, p)[target_name])
    end

    # build optim problem
    optf = Optimization.OptimizationFunction(objective, adtype)
    optprob = Optimization.OptimizationProblem(optf, collect(tunable_pas))
    sol = Optimization.solve(optprob, solve_alg, callback=callback_func, maxiters=maxiters)

    ComponentVector(sol.u, tunable_pas_axes)
end

function param_grad_optimv2(
    component::AbstractComponent;
    tunable_pas::AbstractVector,
    const_pas::AbstractVector,
    input::NamedTuple,
    target::NamedTuple,
    kwargs...,
)
    """
    using Optim.jl for optimization, which has the same problem in param_grad_optim
    """
    # 获取需要优化的参数名称
    target_name = get(kwargs, :target_name, :flow)
    loss_func = get(kwargs, :loss_func, mse)


    tunable_pas_axes = getaxes(tunable_pas)

    # 内部构造一个function
    function predict_func(x::AbstractVector{T}) where {T}
        tmp_tunable_pas = ComponentVector(x, tunable_pas_axes)
        tmp_pas = merge_ca(tmp_tunable_pas, T.(const_pas))[:param]
        output = component(input, tmp_pas)
        println(eltype(output[target_name]))
        output
    end

    objective(x) = loss_func(target[target_name], predict_func(x)[target_name])

    result = Optim.optimize(objective, collect(tunable_pas), LBFGS(); autodiff=:forward)
    result
end

function hybridparams_optimize!()
    """
    混合参数(包括模型超参数和模型内部权重)优化
    """

end

function nn_param_optim(
    nn::AbstractNNFlux;
    input::NamedTuple,
    target::NamedTuple,
    init_params::NamedTuple,
    kwargs...
)
    x = hcat([input[nm] for nm in nn.input_names]...)
    y = hcat([target[nm] for nm in nn.output_names]...)'

    function objective(u, p)
        sum((nn.func(permutedims(x), u)' .- y) .^ 2)
    end

    optf = Optimization.OptimizationFunction(objective, Optimization.AutoZygote())
    optprob = Optimization.OptimizationProblem(optf, ComponentVector(init_params))
    sol = Optimization.solve(optprob, Adam(0.01), maxiters=10)
    sol
end