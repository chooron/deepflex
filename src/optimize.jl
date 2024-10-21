#* callback function for 
default_callback_func(state, l) = begin
    @info (state.iter, l, now())
    false
end

"""
Parameter optimization for global search of hydrological units, nodes
"""
function param_box_optim(
    component::AbstractComponent;
    tunable_pas::ComponentArray,
    const_pas::ComponentArray,
    input::Vector{<:NamedTuple},
    target::Vector{<:NamedTuple},
    timeidx::AbstractVector,
    run_kwargs::Dict=Dict(),
    kwargs...,
)
    #* Get the argument for parameter optimization
    loss_func = get(kwargs, :loss_func, HydroErrors.mse)
    callback_func = get(kwargs, :callback_func, default_callback_func)
    lb = get(kwargs, :lb, zeros(length(tunable_pas)))
    ub = get(kwargs, :ub, ones(length(tunable_pas)) .* 100)
    maxiters = get(kwargs, :maxiters, 10)
    warmup = get(kwargs, :warmup, 100)
    solver = get(kwargs, :solver, ODESolver())
    solve_alg = get(kwargs, :solve_alg, BBO_adaptive_de_rand_1_bin_radiuslimited())

    #* Construct default model parameters based on tunbale parameters and constant parameters for subsequent merging
    default_model_pas = ComponentArray(merge_recursive(NamedTuple(tunable_pas), NamedTuple(const_pas)))
    tunable_axes = getaxes(tunable_pas)

    #* Constructing the objective function for optimization
    function objective(x::AbstractVector{T}, p) where {T}
        #* Optimization arguments: hydro component, input data, time index, ode solver,
        #*                         tunable parameters axes and default model params
        #* Use merge_ca to replace the tunable parameters inner the model parameters
        tmp_pas = update_ca(default_model_pas, ComponentVector(x, tunable_axes))
        loss = mean(map(eachindex(input, target, timeidx)) do i
            tmp_pred = component(input[i], tmp_pas, timeidx[i]; run_kwargs...)
            tmp_loss = mean([loss_func(target[i][key][warmup:end], tmp_pred[key][warmup:end]) for key in keys(target[i])])
            tmp_loss
        end)
        loss
    end

    #* Constructing and solving optimization problems
    optf = Optimization.OptimizationFunction(objective)
    optprob = Optimization.OptimizationProblem(optf, collect(tunable_pas), (), lb=lb, ub=ub)
    sol = Optimization.solve(optprob, solve_alg, callback=callback_func, maxiters=maxiters)
    #* Returns the optimized model parameters
    sol.u
end

function get_objective()
    #* Constructing the objective function for optimization
    function objective(x::AbstractVector{T}, p) where {T}
        #* Optimization arguments: hydro component, input data, time index, ode solver,
        #*                         tunable parameters axes and default model params
        #* Use merge_ca to replace the tunable parameters inner the model parameters
        component, input, target, timeidx, run_kwargs, default_model_pas, loss_func = p
        # todo 添加clip方法,约束参数范围
        tmp_pas = update_ca(default_model_pas, x)
        loss = mean(map(eachindex(input, target, timeidx)) do i
            inp, tar, tidx = input[i], target[i], timeidx[i]
            tmp_pred = component(inp, tmp_pas, tidx; run_kwargs...)
            tmp_loss = mean([loss_func(tar[key], tmp_pred[key]) for key in keys(tar)])
            tmp_loss
        end)
        loss
    end
end

const objective_func = get_objective()

"""
Parameter optimization for local search of hydrological units, nodes
"""
function param_grad_optim(
    component::AbstractComponent;
    tunable_pas::AbstractVector,
    const_pas::AbstractVector,
    input::Vector{<:NamedTuple},
    target::Vector{<:NamedTuple},
    timeidx::AbstractVector,
    run_kwargs::NamedTuple=NamedTuple(),
    kwargs...,
)
    #* Get the argument for parameter optimization
    solve_alg = get(kwargs, :solve_alg, Adam())
    adtype = get(kwargs, :adtype, Optimization.AutoZygote())
    loss_func = get(kwargs, :loss_func, HydroErrors.mse)
    callback_func = get(kwargs, :callback_func, default_callback_func)
    maxiters = get(kwargs, :maxiters, 10)

    #* Construct default model parameters based on tunbale parameters and constant parameters for subsequent merging
    default_model_pas = ComponentArray(merge_recursive(NamedTuple(tunable_pas), NamedTuple(const_pas)))

    #* Constructing and solving optimization problems
    optf = Optimization.OptimizationFunction(objective_func, adtype)
    prob_args = (component, input, target, timeidx, run_kwargs, default_model_pas, loss_func)
    optprob = Optimization.OptimizationProblem(optf, tunable_pas, prob_args)
    sol = Optimization.solve(optprob, solve_alg,  maxiters=maxiters) # callback=callback_func,
    #* Returns the optimized model parameters
    sol.u
end

function batch_param_grad_optim(
    component::AbstractComponent;
    tunable_pas::AbstractVector,
    const_pas::AbstractVector,
    input::Vector{<:AbstractMatrix},
    target::Vector{<:AbstractMatrix},
    timeidx::AbstractVector,
    kwargs...,
)
    #* Get the argument for parameter optimization
    solve_alg = get(kwargs, :solve_alg, Adam())
    adtype = get(kwargs, :adtype, Optimization.AutoZygote())
    loss_func = get(kwargs, :loss_func, HydroErrors.mse)
    maxiters = get(kwargs, :maxiters, 10)
    solver = get(kwargs, :solver, ODESolver())

    train_batch = [(input_i, target_i, timeidx_i) for (input_i, target_i, timeidx_i) in zip(input, target, timeidx)]

    #* Construct default model parameters based on tunbale parameters and constant parameters for subsequent merging
    default_model_pas = ComponentArray(merge_recursive(NamedTuple(tunable_pas), NamedTuple(const_pas)))

    function loss_adjoint(x, _component, _solver, batch_input, batch_target, batch_time)
        tmp_pas = update_ca(default_model_pas, x)
        tmp_pred = _component(batch_input, tmp_pas, timeidx=batch_time, solver=_solver)
        mean([loss_func(batch_target[key], tmp_pred[key]) for key in keys(batch_target)])
    end

    #* Constructing the objective function for optimization
    function objective(x::AbstractVector{T}, p, batch_input, batch_target, batch_time) where {T}
        #* Optimization arguments: hydro component, input data, time index, ode solver,
        #*                         tunable parameters axes and default model params
        #* Use merge_ca to replace the tunable parameters inner the model parameters
        _component, _solver, __default_model_pas = p
        tmp_pas = update_ca(__default_model_pas, x)
        loss = loss_adjoint(tmp_pas, _component, _solver, batch_input, batch_target, batch_time)
        loss
    end

    cumsum_loss = 0.0
    sample_num = length(input)
    batch_callback_func(state, l) = begin
        @info (state.iter, l, now())
        if state.iter % sample_num != 0
            cumsum_loss = cumsum_loss + l
        elseif state.iter % sample_num == 0
            @info Symbol(:epoch_, state.iter ÷ sample_num, Symbol(", mean_loss: "), cumsum_loss / sample_num)
            cumsum_loss = 0.0
        end
        false
    end

    #* Constructing and solving optimization problems
    optf = Optimization.OptimizationFunction(objective, adtype)
    prob_args = (component, solver, default_model_pas)
    optprob = Optimization.OptimizationProblem(optf, tunable_pas, prob_args)
    sol = Optimization.solve(optprob, solve_alg, ncycle(train_batch, maxiters), callback=batch_callback_func)
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
    loss_func = get(kwargs, :loss_func, HydroErrors.mse)
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