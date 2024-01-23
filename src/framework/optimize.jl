function hyper_params_optimize(
    component::C,
    paraminfos::Vector{ParamInfo{T}},
    input::Dict{Symbol,Vector{T}},
    output::Dict{Symbol,Vector{T}},
    weight::Dict{Symbol,T}=Dict{Symbol,T}(),
    errfunc::Dict{Symbol,T}=Dict{Symbol,T}()
) where {C<:Component,T<:Number}
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
    println(sol.u)
    update_paraminfos!(paraminfos, sol.u)
    return paraminfos
end

function hybrid_params_optimize()
    """
    混合参数(包括模型超参数和模型内部权重)优化
    """
end

function nn_parameter_optimize()
    """
    深度学习模型内部参数优化(pretrain)
    """
end