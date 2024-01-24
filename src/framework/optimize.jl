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
    update_paraminfos!(paraminfos, sol.u)
    return paraminfos
end

function hybrid_params_optimize()
    """
    混合参数(包括模型超参数和模型内部权重)优化
    """

end

function nn_params_optimize!(
    ele::LuxElement;
    input::Dict{Symbol,Vector{T}},
    output::Dict{Symbol,Vector{T}},
    epochs::Int=100,
    opt=Adam(0.01f0)) where {T<:Number}
    """
    基于Lux实现深度学习模型内部参数优化(pretrain)
    """
    x = hcat(values(input)...)'
    y = hcat(values(output)...)

    # define loss function
    function loss_function(model, ps, st, data)
        y_pred, st = Lux.apply(model, data[1], ps, st)
        mse_loss = mean(abs2, y_pred .- data[2])
        return mse_loss, st, ()
    end

    # define random seed
    rng = MersenneTwister()
    Random.seed!(rng, 12345)
    tstate = Lux.Training.TrainState(rng, ele.model, opt)
    vjp = Lux.Training.AutoZygote()
    # model training
    (x, y) = (x, y) .|> ele.device
    for epoch in 1:epochs
        grads, loss, stats, tstate = Lux.Training.compute_gradients(vjp,
            loss_function, (x, y), tstate)
        println("Epoch: $(epoch) || Loss: $(loss)")
        tstate = Lux.Training.apply_gradients(tstate, grads)
    end
    update_lux_element!(ele, tstate)
end


function node_params_optimize(
    ele::LuxElement;
    input::Dict{Symbol,Vector{T}},
    output::Dict{Symbol,Vector{T}}) where {T<:Number}
    """
    基于NeuralODE技术实现深度学习模型内部参数优化(pretrain)
    """
    x = ele.device(hcat(values(input)...)')
    y = ele.device(hcat(values(output)...))

    function loss_function(ps, x, y)
        pred, st_ = ele.model(x, ps, ele.states)
        return mse(pred, y), pred
    end
    opt_func = OptimizationFunction((ps, _, x, y) -> loss_function(ps, x, y), Optimization.AutoZygote())
    opt_prob = OptimizationProblem(opt_func, ele.parameters)
    res = Optimization.solve(opt_prob, opt, zip(x_train, y_train); callback)
end