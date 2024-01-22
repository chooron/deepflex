struct Optimize
    optmzr

end

function mse(predict, target)
    sum(abs(target .- predict)) / length(target)
end

function hyper_params_optimize(
    component::C,
    paraminfos::Vector{ParamInfo},
    input::Dict{Symbol,Vector{T}},
    output::Dict{Symbol,Vector{T}},
    weight::Dict{Symbol,T}=Dict(),
    errfunc::Dict{Symbol,T}=Dict()
) where {C<:Component,T<:Number}
    """
    针对模型超参数进行优化
    """
    # 设置默认weight和errfunc
    if length(weight) == 0
        weight = Dict(1.0 for k in keys(target))
    end
    if length(weight) == 0
        errfunc = Dict(mse for k in keys(target))
    end

    # 内部构造一个function
    function objective(x, p)
        # 为所有componet创建set_parameters函数
        update_paraminfos!(paraminfos, x)
        set_parameters!(component, paraminfos=paraminfos)
        predict = get_output(component, input=input)
        criteria = 0.0
        for (k, v) in output
            criteria += err_func[k](v, predict[k]) * weight[k]
        end
        return criteria
    end
    optf = Optimization.OptimizationFunction(objective, Optimization.AutoZygote())
    optprob = Optimization.OptimizationProblem(optf, p_init)
    sol = Optimization.solve(optprob, optmzr, callback=callback, maxiters=max_N_iter)
    refine_paraminfos!(paraminfos, sol.u)
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