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

struct HyperOptimizer <: AbstractOptimizer
    optf
    optprob
    alg
end

function HyperOptimizer(; func::Function, paraminfos::Vector{BPI}, kwargs...) where {BPI<:BoundaryParamInfo}
    x0 = [p.default for p in paraminfos]
    lb = [p.lb for p in paraminfos]
    ub = [p.ub for p in paraminfos]
    optf = Optimization.OptimizationFunction(func)
    optprob = Optimization.OptimizationProblem(optf, x0, (), lb=lb, ub=ub)

    HyperOptimizer(
        optf,
        optprob,
        kwargs[:alg]
    )
end

function opt_solve(optimizer::HyperOptimizer; kwargs...)
    if haskey(kwargs, :callback)
        callback = kwargs[:callback]
    else
        callback = function (p, l)
            println("rmse: " * string(l))
            return false
        end
    end
    Optimization.solve(optimizer.optprob, optimizer.alg, callback=callback, maxiters=kwargs[:maxiters])
end

function hyperparams_optimize(
    component::AbstractComponent;
    paraminfos::Vector{P},
    input::ComponentVector{T},
    target::ComponentVector{T},
    kwargs...,
) where {P<:AbstractParamInfo,T<:Number}
    """
    针对模型超参数进行优化
    """
    # 设置默认weight和errfunc
    if haskey(kwargs, :error_weight)
        error_weight = kwargs[:error_weight]
    else
        error_weight = Dict(k => 1.0 for k in keys(target))
    end

    if haskey(kwargs, :error_func)
        error_func = kwargs[:error_func]
    else
        error_func = Dict(k => rmse for k in keys(target))
    end

    param_names = [p.name for p in paraminfos]

    # 内部构造一个function
    function objective(x, p)
        tmp_params = ComponentVector(namedtuple(param_names, x))
        predict = get_output(component, input=input, parameters=tmp_params,
            init_states=tmp_params[component.state_names])
        criteria = 0.0
        for k in keys(target)
            criteria += error_func[k](target[k], predict[k]) * error_weight[k]
        end
        return criteria
    end

    hyper_optimizer = HyperOptimizer(func=objective, paraminfos=paraminfos, alg=BBO_adaptive_de_rand_1_bin_radiuslimited())
    opt_solve(hyper_optimizer, maxiters=100)
end

function hybridparams_optimize!()
    """
    混合参数(包括模型超参数和模型内部权重)优化
    """

end

function pretrain!(nn::LuxNNFlux; input::ComponentVector{T}, train_config...) where {T<:Number}
    x = hcat([input[nm] for nm in nn.input_names]...)
    y = hcat([input[nm] for nm in nn.output_names]...)'

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
    nn.parameters = ComponentArray(nn.parameters; Dict(:ps => sol.u)...)
end
