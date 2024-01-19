struct Optimize
    optmzr

end

function hyper_params_optimize(component::C, parameters <: Dict{Symbol,T}; optmzr) where {C<:Component,T<:Number}
    """
    针对模型超参数进行优化
    """
    # 内部构造一个function
    function objective(p)

    end
    optf = Optimization.OptimizationFunction((θ, p) -> loss_model(θ), Optimization.AutoZygote())
    optprob = Optimization.OptimizationProblem(optf, p_init)
    sol = Optimization.solve(optprob, optmzr, callback = callback, maxiters = max_N_iter)

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