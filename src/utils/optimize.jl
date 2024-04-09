function param_box_optim(
    component::AbstractComponent;
    search_params::AbstractVector,
    const_params::AbstractVector,
    input::NamedTuple,
    target::NamedTuple,
    kwargs...,
)
    """
    针对模型参数进行二次优化
    """
    # 获取需要优化的参数名称
    search_param_names = [p.name for p in search_params]
    const_param_names = [p.name for p in const_params]
    const_param_values = [p.value for p in const_params]
    default_values = [p.default for p in search_params]
    lb = [p.lb for p in search_params]
    ub = [p.ub for p in search_params]

    # 设置默认weight和errfunc
    error_weight = get(kwargs, :error_weight, Dict(k => 1.0 for k in keys(target)))
    error_func = get(kwargs, :error_func, Dict(k => rmse for k in keys(target)))
    solve_alg = get(kwargs, :solve_alg, BBO_adaptive_de_rand_1_bin_radiuslimited())

    # TODO 联合callback库有一个统一的标准
    callback_func = get(kwargs, :callback_func, (p, l) -> begin
        @info l
        false
    end)

    # 内部构造一个function
    function objective(x::AbstractVector{T}, p) where {T}
        tmp_search_params = namedtuple(search_param_names, x)
        search_params_type = eltype(tmp_search_params)
        tmp_params = merge(tmp_search_params, namedtuple(const_param_names, search_params_type.(const_param_values)))
        predict = component(input, tmp_params[component.nameinfo.param_names], tmp_params[component.nameinfo.state_names])
        sum((target[:flow] .- predict[:flow]) .^ 2) / length(target[:flow])
    end

    # 构建问题
    optf = Optimization.OptimizationFunction(objective)
    optprob = Optimization.OptimizationProblem(optf, default_values, lb=lb, ub=ub)
    sol = Optimization.solve(optprob, solve_alg, callback=callback_func, maxiters=get(kwargs, :maxiters, 10))
    namedtuple(search_param_names, sol.u)
end

function param_grad_optim(
    component::AbstractComponent;
    search_params::AbstractVector,
    const_params::AbstractVector,
    input::NamedTuple,
    target::NamedTuple,
    kwargs...,
)
    """
    针对模型参数进行二次优化
    """
    # 获取需要优化的参数名称
    search_param_names = [p.name for p in search_params]
    default_values = [p.default for p in search_params]

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
        tmp_search_params = namedtuple(search_param_names, u)
        tmp_params = merge(tmp_search_params, namedtuple([p.name for p in const_params], eltype(u).([p.value for p in const_params])))
        predict = component(input, tmp_params[component.nameinfo.param_names], tmp_params[component.nameinfo.state_names])
        loss = sum((target[:flow] .- predict[:flow]) .^ 2) / length(target[:flow])
        return loss
    end

    # 构建问题
    optf = Optimization.OptimizationFunction(objective, get(kwargs, :adtype, Optimization.AutoForwardDiff())) # AutoForwardDiff and AutoFiniteDiff
    optprob = Optimization.OptimizationProblem(optf, default_values)
    sol = Optimization.solve(optprob, solve_alg, callback=callback_func, maxiters=maxiters)
    namedtuple(search_param_names, sol.u)
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