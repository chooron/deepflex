#* callback function for 
default_callback_func(p, l) = begin
    @info l
    false
end

#* build predict function
function predict_func(x::AbstractVector{T}, p) where {T}
    #* Optimization arguments: hydro component, input data, time index, ode solver,
    #*                         tunable parameters axes and default model params
    component, input, timeidx, solver, default_model_pas = p
    #* Use merge_ca to replace the tunable parameters inner the model parameters
    tmp_pas = merge_ca(default_model_pas, x)
    #* Call the hydro component to calculate the simulation results
    component(input, tmp_pas, timeidx=timeidx, solver=solver)
end

"""
$(SIGNATURES)

Parameter optimization for global search of hydrological units, nodes
"""
function param_box_optim(
    component::AbstractComponent;
    tunable_pas::ComponentArray,
    const_pas::ComponentArray,
    input::NamedTuple,
    target::NamedTuple,
    timeidx::Vector,
    kwargs...,
)
    #* Get the argument for parameter optimization
    loss_func = get(kwargs, :loss_func, mse)
    callback_func = get(kwargs, :callback_func, default_callback_func)
    lb = get(kwargs, :lb, zeros(length(tunable_pas)))
    ub = get(kwargs, :ub, ones(length(tunable_pas)) .* 100)
    maxiters = get(kwargs, :maxiters, 10)
    solver = get(kwargs, :solver, ODESolver(saveat=timeidx))
    solve_alg = get(kwargs, :solve_alg, BBO_adaptive_de_rand_1_bin_radiuslimited())

    #* Construct default model parameters based on tunbale parameters and constant parameters for subsequent merging
    default_model_pas = ComponentArray(merge_recursive(NamedTuple(tunable_pas), NamedTuple(const_pas)))
    #* Get the axes of the tunbale parameters
    tunable_pas_axes = getaxes(tunable_pas)

    #* Constructing the objective function for optimization
    function objective(x::AbstractVector{T}, p) where {T}
        predict_result = predict_func(x, p)
        loss =  mean([loss_func(target[key], predict_result[key]) for key in keys(target)])
        loss
    end

    #* Constructing and solving optimization problems
    optf = Optimization.OptimizationFunction(objective)
    prob_args = (component, input, timeidx, solver, tunable_pas_axes, default_model_pas)
    optprob = Optimization.OptimizationProblem(optf, collect(tunable_pas), prob_args, lb=lb, ub=ub)
    sol = Optimization.solve(optprob, solve_alg, callback=callback_func, maxiters=maxiters)
    #* Returns the optimized model parameters
    ComponentVector(sol.u, tunable_pas_axes)
end

"""
$(SIGNATURES)

Parameter optimization for local search of hydrological units, nodes
"""
function param_grad_optim(
    component::AbstractComponent;
    tunable_pas::AbstractVector,
    const_pas::AbstractVector,
    input::NamedTuple,
    target::NamedTuple,
    timeidx::Vector,
    kwargs...,
)
    #* Get the argument for parameter optimization
    solve_alg = get(kwargs, :solve_alg, Adam())
    adtype = get(kwargs, :adtype, Optimization.AutoZygote())
    loss_func = get(kwargs, :loss_func, mse)
    callback_func = get(kwargs, :callback_func, default_callback_func)
    maxiters = get(kwargs, :maxiters, 10)
    solver = get(kwargs, :solver, ODESolver())

    #* Construct default model parameters based on tunbale parameters and constant parameters for subsequent merging
    default_model_pas = ComponentArray(merge_recursive(NamedTuple(tunable_pas), NamedTuple(const_pas)))
    # #* Get the axes of the tunbale parameters
    # tunable_pas_axes = getaxes(tunable_pas)

    #* Constructing the objective function for optimization
    function objective(x::AbstractVector{T}, p) where {T}
        predict_result = predict_func(x, p)
        loss = mean([loss_func(target[key], predict_result[key]) for key in keys(target)])
        loss
    end

    #* Constructing and solving optimization problems
    optf = Optimization.OptimizationFunction(objective, adtype)
    prob_args = (component, input, timeidx, solver, default_model_pas)
    optprob = Optimization.OptimizationProblem(optf, tunable_pas, prob_args)
    sol = Optimization.solve(optprob, solve_alg, callback=callback_func, maxiters=maxiters)
    #* Returns the optimized model parameters
    sol.u
end

function nn_param_optim(
    flux::AbstractNeuralFlux;
    input::NamedTuple,
    target::NamedTuple,
    init_params::ComponentVector,
    kwargs...
)
    #* Get the argument for parameter optimization
    solve_alg = get(kwargs, :solve_alg, Adam(0.01))
    adtype = get(kwargs, :adtype, Optimization.AutoZygote())
    maxiters = get(kwargs, :maxiters, 100)
    loss_func = get(kwargs, :loss_func, mse)
    callback_func = get(kwargs, :callback_func, default_callback_func)

    #* Integrate nn's input variables
    x = reduce(hcat, collect(input[get_input_names(flux)]))
    flux_nn_name = get_param_names(flux)[1]
    flux_output_name = get_output_names(flux)

    #* Constructing the objective function for optimization
    function objective(u, p)
        predict = flux(x, [u[flux_nn_name]])
        predict_ntp = NamedTuple{Tuple(flux_output_name)}(eachcol(predict))
        mean([loss_func(predict_ntp[nm], target[nm]) for nm in keys(target)])
    end

    #* Constructing and solving optimization problems
    optf = Optimization.OptimizationFunction(objective, adtype)
    optprob = Optimization.OptimizationProblem(optf, init_params)
    sol = Optimization.solve(optprob, solve_alg,callback=callback_func, maxiters=maxiters)
    sol.u
end