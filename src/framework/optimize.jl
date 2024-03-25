using ForwardDiff

mutable struct BoundaryParamInfo{T} <: AbstractParamInfo where {T<:Number}
    name::Symbol
    default::T
    lb::T
    ub::T
    value::T
end

function BoundaryParamInfo(name::Symbol, default::T;
    lb::T=default, ub::T=default, value::T=default) where {T<:Number}
    BoundaryParamInfo(
        name,
        default,
        lb,
        ub,
        value
    )
end

function hyperparams_gradient_optimize(
    component::AbstractComponent;
    search_params::Vector{P},
    const_params::NamedTuple,
    input::NamedTuple,
    target::NamedTuple,
    kwargs...,
) where {P<:AbstractParamInfo}
    """
    针对模型超参数进行优化
    """
    # 获取需要优化的参数名称
    search_param_names = [p.name for p in search_params]

    # 设置默认weight和errfunc
    error_weight = get(kwargs, :error_weight, Dict(k => 1.0 for k in keys(target)))
    error_func = get(kwargs, :error_func, Dict(k => rmse for k in keys(target)))
    solve_alg = get(kwargs, :solve_alg, Adam())
    callback_func = get(kwargs, :callback_func, (p, l) -> begin
        println("loss: " * string(l))
        return false
    end)

    # 内部构造一个function
    function objective(x::AbstractVector{T}, p)  where T
        tmp_search_params = namedtuple(search_param_names, x)
        search_params_type = eltype(tmp_search_params)
        tmp_params = merge(tmp_search_params, namedtuple(keys(const_params), search_params_type.(collect(const_params))))
        predict = get_output(component, input=input, parameters=tmp_params,
            init_states=tmp_params[component.state_names])
        sum((target[:flow] .- predict[:flow]) .^ 2) / length(target[:flow])
    end

    # 处理数据搜索范围以及默认值
    x0 = [p.default for p in search_params]
    lb = [p.lb for p in search_params]
    ub = [p.ub for p in search_params]
    # 构建问题
    optf = Optimization.OptimizationFunction(objective, get(kwargs, :adtype, Optimization.AutoForwardDiff())) # get(kwargs, :adtype, Optimization.AutoZygote())
    optprob = Optimization.OptimizationProblem(optf, x0) # , lb=lb, ub=ub
    sol = Optimization.solve(optprob, solve_alg, callback=callback_func, maxiters=get(kwargs, :maxiters, 10))
    namedtuple(search_param_names, sol.u)
end

function hybridparams_optimize!()
    """
    混合参数(包括模型超参数和模型内部权重)优化
    """

end

function nn_params_optimize!(nn::LuxNNFlux; input::NamedTuple, output::NamedTuple, train_config...)
    x = hcat([input[nm] for nm in nn.input_names]...)
    y = hcat([output[nm] for nm in nn.output_names]...)'

    function prep_pred_NN_pretrain(model_, input_)
        (params) -> model_(input_, params)
    end

    pred_NN_pretrain_fct = prep_pred_NN_pretrain(nn.func, permutedims(x))

    function loss_NN_pretrain(params, batch)
        sum((pred_NN_pretrain_fct(params)' .- batch) .^ 2)
    end

    optf = Optimization.OptimizationFunction((θ, p) -> loss_NN_pretrain(θ, y), Optimization.AutoZygote())
    optprob = Optimization.OptimizationProblem(optf, nn.parameters[:ps])
    sol = Optimization.solve(optprob, Adam(0.01), maxiters=100)
    ComponentArray(nn.parameters; Dict(:ps => sol.u)...)
end
