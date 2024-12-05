
function param_batch_optim(
    component::AbstractComponent;
    tunable_pas::ComponentVector,
    const_pas::ComponentVector,
    input::Vector,
    target::Vector,
    config::Vector{<:NamedTuple}=[NamedTuple()],
    run_kwargs::NamedTuple=(convert_to_ntp=true,),
    opt_kwargs...,
)
    #* Get the argument for parameter optimization
    loss_func = get(opt_kwargs, :loss_func, (obs, sim) -> sum((obs .- sim) .^ 2) / length(obs))
    adtype = get(opt_kwargs, :adtype, Optimization.AutoZygote())
    maxiters = get(opt_kwargs, :maxiters, 100)
    warmup = get(opt_kwargs, :warmup, 100)
    solve_alg = get(opt_kwargs, :solve_alg, Adam())

    #* prepare the batch data
    train_batch = [(input_i, target_i, cfg_i) for (input_i, target_i, cfg_i) in zip(input, target, config)]

    #* Construct default model parameters based on tunbale parameters and constant parameters for subsequent merging
    default_model_pas = ComponentArray(merge_recursive(NamedTuple(tunable_pas), NamedTuple(const_pas)))

    #* Constructing the objective function for optimization
    function objective(x::AbstractVector{T}, p, batch_input, batch_target, batch_cfg) where {T}
        #* Optimization arguments: hydro component, input data, time index, ode solver,
        #*                         tunable parameters axes and default model params
        #* Use merge_ca to replace the tunable parameters inner the model parameters
        _component, __default_model_pas = p
        tmp_pas = update_ca(__default_model_pas, x)
        tmp_pred = _component(batch_input, tmp_pas; config=batch_cfg, run_kwargs...)
        loss = mean([loss_func(batch_target[key][warmup:end], tmp_pred[key][warmup:end]) for key in keys(batch_target)])
        loss
    end

    # todo 添加earlystop功能
    # if all([k in keys(opt_kwargs) for k in [:val_input, :val_target, :val_timeidx]])
    #     #* validation arguments
    #     patience = get(opt_kwargs, :patience, 10)
    #     min_epoch = get(opt_kwargs, :min_epoch, 5)
    #     val_input = get(opt_kwargs, :val_input, input)
    #     val_target = get(opt_kwargs, :val_target, target)
    #     val_timeidx = get(opt_kwargs, :val_timeidx, timeidx)
    #     val_dataset = [val_input, val_target, val_timeidx]
    #     loss_recorder = NamedTuple[]
    #     callback_func = get(opt_kwargs, :callback_func, get_batch_callback_func_with_earlystop(length(input), loss_recorder, val_dataset, patience, min_epoch))
    loss_recorder = NamedTuple[]
    callback_func = get(opt_kwargs, :callback_func, get_batch_callback_func(length(input), loss_recorder))

    #* Constructing and solving optimization problems
    optf = Optimization.OptimizationFunction(objective, adtype)
    optprob = Optimization.OptimizationProblem(optf, tunable_pas, (component, default_model_pas))
    sol = Optimization.solve(optprob, solve_alg, ncycle(train_batch, maxiters), callback=callback_func)
    #* Returns the optimized model parameters
    sol.u
end

function nn_param_optim(
    flux::AbstractNeuralFlux;
    input::AbstractMatrix,
    target::NamedTuple,
    init_params::ComponentVector,
    kwargs...
)
    #* Get the argument for parameter optimization
    solve_alg = get(kwargs, :solve_alg, Adam(0.01))
    adtype = get(kwargs, :adtype, Optimization.AutoZygote())
    maxiters = get(kwargs, :maxiters, 100)
    loss_func = get(kwargs, :loss_func, (obs, sim) -> sum((obs .- sim) .^ 2) / length(obs))
    callback_func = get(kwargs, :callback_func, default_callback_func)

    #* Integrate nn's input variables
    flux_nn_name = get_param_names(flux)[1]
    flux_output_name = get_output_names(flux)

    #* Constructing the objective function for optimization
    function objective(u, p)
        predict = flux(input, [u[flux_nn_name]])
        predict_ntp = NamedTuple{Tuple(flux_output_name)}(eachcol(predict))
        mean([loss_func(predict_ntp[nm], target[nm]) for nm in keys(target)])
    end

    #* Constructing and solving optimization problems
    optf = Optimization.OptimizationFunction(objective, adtype)
    optprob = Optimization.OptimizationProblem(optf, init_params)
    sol = Optimization.solve(optprob, solve_alg, callback=callback_func, maxiters=maxiters)
    sol.u
end