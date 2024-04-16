mse(y, y_hat) = sum((y .- y_hat) .^ 2) / length(y)
default_callback_func(p, l) = begin
    @info l
    false
end

function param_box_optim(
    component::AbstractComponent;
    tunable_params::ComponentArray,
    const_params::ComponentArray,
    input::NamedTuple,
    target::NamedTuple,
    kwargs...,
)
    """
    针对模型参数进行二次优化
    """
    # 获取需要优化的参数名称
    tunable_param_axes = getaxes(tunable_params)

    solve_alg = get(kwargs, :solve_alg, BBO_adaptive_de_rand_1_bin_radiuslimited())
    target_name = get(kwargs, :target_name, :flow)
    loss_func = get(kwargs, :loss_func, mse)
    callback_func = get(kwargs, :callback_func, default_callback_func)
    lb = get(kwargs, :lb, zeros(length(tunable_params)))
    ub = get(kwargs, :ub, ones(length(tunable_params)) .* 100)
    maxiters=get(kwargs, :maxiters, 10)

    predict_func(x, p) = begin
        tunable_params_type = eltype(x)
        tmp_tunable_params = ComponentVector(x, tunable_param_axes)
        tmp_params = ComponentVector(tmp_tunable_params; tunable_params_type.(const_params)...)
        component(input, tmp_params, tmp_params)
    end

    objective(x, p) = loss_func(target[target_name], predict_func(x, p)[target_name])

    # 构建问题
    optf = Optimization.OptimizationFunction(objective)
    optprob = Optimization.OptimizationProblem(optf, collect(tunable_params), lb=lb, ub=ub)
    sol = Optimization.solve(optprob, solve_alg, callback=callback_func, maxiters=maxiters)
    namedtuple(tunable_param_names, sol.u)
end

function param_grad_optim(
    component::AbstractComponent;
    tunable_params::AbstractVector,
    const_params::AbstractVector,
    input::NamedTuple,
    target::NamedTuple,
    kwargs...,
)
    """
    针对模型参数进行二次优化
    """
    # 获取需要优化的参数名称
    tunable_param_names = [p.name for p in tunable_params]
    default_values = [p.default for p in tunable_params]

    # 设置默认weight和errfunc
    error_weight = get(kwargs, :error_weight, Dict(k => 1.0 for k in keys(target)))
    error_func = get(kwargs, :error_func, Dict(k => rmse for k in keys(target)))

    solve_alg = get(kwargs, :solve_alg, Adam())
    maxiters = get(kwargs, :maxiters, 10)

    # build problem
    # for ele in vcat(component.surf_layer, component.soil_layer, component.slope_layer)
    #     if length(ele.dfuncs) > 0
    #         build_prob!(ele, input=input)
    #     end
    # end

    # TODO 联合callback库有一个统一的标准
    function callback_func(p, l)
        @info l
        false
    end

    # 内部构造一个function
    function objective(u, p)
        tmp_tunable_params = namedtuple(tunable_param_names, u)
        tmp_params = merge(tmp_tunable_params, namedtuple([p.name for p in const_params], eltype(u).([p.value for p in const_params])))
        predict = component(input, tmp_params[component.param_names], tmp_params[component.state_names])
        loss = sum((target[:flow] .- predict[:flow]) .^ 2) / length(target[:flow])
        return loss
    end

    # 构建问题
    optf = Optimization.OptimizationFunction(objective, get(kwargs, :adtype, Optimization.AutoForwardDiff())) # AutoForwardDiff and AutoFiniteDiff
    optprob = Optimization.OptimizationProblem(optf, default_values)
    sol = Optimization.solve(optprob, solve_alg, callback=callback_func, maxiters=maxiters)
    namedtuple(tunable_param_names, sol.u)
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