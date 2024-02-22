function hyper_params_optimize(
    component::C,
    paraminfos::Vector{P},
    input::Dict{Symbol,Vector{T}},
    output::Dict{Symbol,Vector{T}},
    weight::Dict{Symbol,T}=Dict{Symbol,T}(),
    errfunc::Dict{Symbol,T}=Dict{Symbol,T}()
) where {C<:AbstractComponent,P<:AbstractParamInfo,T<:Number}
    """
    针对模型超参数进行优化
    """
    # 设置默认weight和errfunc
    if length(weight) == 0
        weight = Dict(k => 1.0 for k in keys(output))
    end
    if length(errfunc) == 0
        errfunc = Dict(k => rmse for k in keys(output))
    end

    # 内部构造一个function
    function objective(x, p)
        # 为所有componet创建set_parameters函数
        update_paraminfos!(paraminfos, x)
        set_parameters!(component, paraminfos=paraminfos)
        predict = get_output(component, input=input)
        criteria = 0.0
        for (k, v) in output
            criteria += errfunc[k](v, predict[k]) * weight[k]
        end
        return criteria
    end

    callback = function (p, l)
        println("rmse: " * string(l))
        return false
    end

    x0 = [p.default for p in paraminfos]
    lb = [p.lb for p in paraminfos]
    ub = [p.ub for p in paraminfos]
    optf = Optimization.OptimizationFunction(objective)
    optprob = Optimization.OptimizationProblem(optf, x0, (), lb=lb, ub=ub)
    sol = Optimization.solve(optprob, BBO_adaptive_de_rand_1_bin_radiuslimited(), callback=callback, maxiters=100)
    update_paraminfos!(paraminfos, sol.u)
    return paraminfos
end

function hybrid_params_optimize()
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
