


function get_objective()
    #* Constructing the objective function for optimization
    function objective(x::AbstractVector{T}, p) where {T}
        #* Optimization arguments: hydro component, input data, time index, ode solver,
        #*                         tunable parameters axes and default model params
        component, input, target, config, run_kwargs, default_model_pas, loss_func, warmup = p
        #* Use merge_ca to replace the tunable parameters inner the model parameters
        tmp_pas = update_ca(default_model_pas, x)
        loss = mean(map(eachindex(input, target, config)) do i
            inp, tar, cfg = input[i], target[i], config[i]
            tmp_pred = component(inp, tmp_pas; config=cfg, run_kwargs...)
            tmp_loss = mean([loss_func(tar[key][warmup:end], tmp_pred[key][warmup:end]) for key in keys(tar)])
            tmp_loss
        end)
        loss
    end

    return objective
end

const objective_func = get_objective()

"""
    param_box_optim(component::AbstractComponent; kwargs...)

Perform box-constrained parameter optimization for hydrological units or nodes.

# Arguments
- `component::AbstractComponent`: The hydrological component to optimize.
- `tunable_pas::ComponentVector`: Vector of tunable parameters.
- `const_pas::ComponentVector`: Vector of constant parameters.
- `input::Vector{<:NamedTuple}`: Vector of input data for each optimization iteration.
- `target::Vector{<:NamedTuple}`: Vector of target data for each optimization iteration.
- `timeidx::AbstractVector`: Vector of time indices for each optimization iteration.
- `run_kwargs::Dict=Dict()`: Additional keyword arguments for running the component.

# Optional keyword arguments
- `loss_func::Function=HydroErrors.mse`: Loss function to use for optimization.
- `maxiters::Int=10`: Maximum number of iterations for optimization.
- `warmup::Int=100`: Number of warmup steps before calculating loss.
- `solve_alg::Optimization.AbstractOptimizer=BBO_adaptive_de_rand_1_bin_radiuslimited()`: Optimization algorithm to use.
- `lb::Vector`: Lower bounds for parameters (default: zeros).
- `ub::Vector`: Upper bounds for parameters (default: 100 for each parameter).
- `callback_func::Function`: Custom callback function for optimization process.

# Returns
- `sol.u`: Optimized parameters.
- `loss_df`: A DataFrame containing the loss history.

# Description
This function performs box-constrained parameter optimization for a given hydrological component. 
It uses a black-box optimization algorithm to find the best set of parameters within specified bounds 
that minimize the loss between the model output and the target data.
"""
function param_box_optim(
    component::AbstractComponent;
    tunable_pas::ComponentVector,
    const_pas::ComponentVector,
    input::Vector{<:NamedTuple},
    target::Vector{<:NamedTuple},
    config::Vector{<:NamedTuple}=fill(NamedTuple(), length(input)),
    run_kwargs::NamedTuple=(convert_to_ntp=true,),
    opt_kwargs...,
)
    @assert length(input) == length(target) == length(config) "The length of input, target and config must be the same,
     while $(length(input)) input, $(length(target)) target, $(length(config)) config are given."
    #* Get the argument for parameter optimization
    loss_func = get(opt_kwargs, :loss_func, (obs, sim) -> sum((obs .- sim) .^ 2) / length(obs))
    loss_recorder = NamedTuple[]
    lb = get(opt_kwargs, :lb, zeros(length(tunable_pas)))
    ub = get(opt_kwargs, :ub, ones(length(tunable_pas)) .* 100)
    maxiters = get(opt_kwargs, :maxiters, 10)
    warmup = get(opt_kwargs, :warmup, 100)
    callback_func = get(opt_kwargs, :callback_func, get_callback_func(Progress(maxiters, desc="Training..."), loss_recorder))
    solve_alg = get(opt_kwargs, :solve_alg, BBO_adaptive_de_rand_1_bin_radiuslimited())
    tunable_axes = getaxes(tunable_pas)
    #* Construct default model parameters based on tunbale parameters and constant parameters for subsequent merging
    default_model_pas = ComponentArray(merge_recursive(NamedTuple(tunable_pas), NamedTuple(const_pas)))

    #* Constructing the objective function for optimization
    function objective(x::AbstractVector{T}, p) where {T}
        tmp_pas = update_ca(default_model_pas, ComponentVector(x, tunable_axes))
        loss = mean(map(eachindex(input, target, config)) do i
            inp, tar, cfg = input[i], target[i], config[i]
            tmp_pred = component(inp, tmp_pas; config=cfg, run_kwargs...)
            tmp_loss = mean([loss_func(tar[key][warmup:end], tmp_pred[key][warmup:end]) for key in keys(tar)])
            tmp_loss
        end)
        loss
    end

    #* Constructing and solving optimization problems
    optf = Optimization.OptimizationFunction(objective)
    optprob = Optimization.OptimizationProblem(optf, collect(tunable_pas), (), lb=lb, ub=ub)
    sol = Optimization.solve(optprob, solve_alg, callback=callback_func, maxiters=maxiters)
    #* convert loss_recorder to DataFrame
    loss_df = DataFrame(loss_recorder)
    #* Returns the optimized model parameters and the loss recorder
    sol.u, loss_df
end

"""
    param_grad_optim(component::AbstractComponent; kwargs...)

Perform gradient-based parameter optimization for local search of hydrological units or nodes.

# Arguments
- `component::AbstractComponent`: The hydrological component to optimize.
- `tunable_pas::ComponentVector`: Vector of tunable parameters.
- `const_pas::ComponentVector`: Vector of constant parameters.
- `input::Vector{<:NamedTuple}`: Vector of input data for each optimization iteration.
- `target::Vector{<:NamedTuple}`: Vector of target data for each optimization iteration.
- `timeidx::AbstractVector`: Vector of time indices for each optimization iteration.
- `run_kwargs::NamedTuple=NamedTuple()`: Additional keyword arguments for running the component.

# Optional keyword arguments
- `loss_func::Function=HydroErrors.mse`: Loss function to use for optimization.
- `adtype::Optimization.AbstractADType=Optimization.AutoZygote()`: Automatic differentiation type.
- `maxiters::Int=10`: Maximum number of iterations for optimization.
- `warmup::Int=100`: Number of warmup steps before calculating loss.
- `solve_alg::Optimization.AbstractOptimizer=Adam()`: Optimization algorithm to use.
- `callback_func::Function`: Custom callback function for optimization process.

# Returns
- `sol.u`: Optimized parameters.
- `loss_df`: A DataFrame containing the loss history.

# Description
This function performs gradient-based parameter optimization for a given hydrological component. 
It uses automatic differentiation to compute gradients and applies an optimization algorithm 
to find the best set of parameters that minimize the loss between the model output and the target data.
"""
function param_grad_optim(
    component::AbstractComponent;
    tunable_pas::ComponentVector,
    const_pas::ComponentVector,
    input::Vector{<:NamedTuple},
    target::Vector{<:NamedTuple},
    config::Vector{<:NamedTuple}=fill(NamedTuple(), length(input)),
    run_kwargs::NamedTuple=(convert_to_ntp=true,),
    opt_kwargs...,
)
    @assert(length(input) == length(target) == length(config),
        "The length of input, target and config must be the same, while $(length(input)) input, $(length(target)) target, $(length(config)) config are given.")

    #* Get the argument for parameter optimization
    loss_func = get(opt_kwargs, :loss_func, (obs, sim) -> sum((obs .- sim) .^ 2) / length(obs))
    adtype = get(opt_kwargs, :adtype, Optimization.AutoZygote())
    maxiters = get(opt_kwargs, :maxiters, 10)
    warmup = get(opt_kwargs, :warmup, 100)
    solve_alg = get(opt_kwargs, :solve_alg, Adam())

    loss_recorder = NamedTuple[]
    callback_func = get(opt_kwargs, :callback_func, get_callback_func(Progress(maxiters, desc="Training..."), loss_recorder))

    #* Construct default model parameters based on tunbale parameters and constant parameters for subsequent merging
    default_model_pas = ComponentArray(merge_recursive(NamedTuple(tunable_pas), NamedTuple(const_pas)))

    #* Constructing and solving optimization problems
    optf = Optimization.OptimizationFunction(objective_func, adtype)
    prob_args = (component, input, target, config, run_kwargs, default_model_pas, loss_func, warmup)
    optprob = Optimization.OptimizationProblem(optf, tunable_pas, prob_args)
    sol = Optimization.solve(optprob, solve_alg, callback=callback_func, maxiters=maxiters)
    #* convert loss_recorder to DataFrame
    loss_df = DataFrame(loss_recorder)
    #* Returns the optimized model parameters and the loss recorder
    sol.u, loss_df
end

function param_batch_optim(
    component::AbstractComponent;
    tunable_pas::ComponentVector,
    const_pas::ComponentVector,
    input::Vector{<:NamedTuple},
    target::Vector{<:NamedTuple},
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