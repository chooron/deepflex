"""
$(TYPEDEF)
The basic hydrological calculation module contains multiple hydrological fluxes,
and can simulate the balance calculation of a physical module.
# Fields
$(FIELDS)
# Example
```
```
"""
struct HydroElement <: AbstractHydroElement
    "the name of hydrological computation element"
    name::Symbol
    "element information: keys contains: input, output, param, state"
    nameinfo::NamedTuple
    "common hydrological fluxes, used to provide calculation results for state fluxes"
    funcs::Vector
    """
    Hydrological state fluxes, 
    combined with ordinary hydrological flux to construct ordinary differential equations
    """
    dfuncs::Vector
    """
    Hydrological flux functions
    """
    flux_func::Function
    """
    Hydrological ode functions
    """
    ode_func::Union{Nothing,Function}

    function HydroElement(
        name::Symbol;
        funcs::Vector,
        dfuncs::Vector=StateFlux[],
    )
        #* Extract all variable names of funcs and dfuncs
        ele_input_names, ele_output_names, ele_state_names = get_var_names(funcs, dfuncs)
        #* Extract all parameters names of funcs and dfuncs
        ele_param_names = get_param_names(vcat(funcs, dfuncs))
        #* Extract all neuralnetwork names of the funcs
        ele_nn_names = get_nn_names(funcs)
        #* Setup the name information of the hydroelement
        nameinfo = (input=ele_input_names, output=ele_output_names, state=ele_state_names, param=ele_param_names, nn=ele_nn_names)
        #* Construct a NamedTuple of all func's variables and parameters
        funcs_vars_ntp = reduce(merge, [merge(func.input_info, func.output_info) for func in vcat(funcs, dfuncs)])
        # todo 这块尽量不要包括
        funcs_params_ntp = reduce(merge, [func.param_info for func in vcat(funcs, dfuncs) if !(func isa AbstractNeuralFlux)])
        funcs_nn_params = collect([func.nn_info[:params] for func in funcs if func isa AbstractNeuralFlux])
        #* Construct a function for ordinary differential calculation based on dfunc and funcs
        flux_func, ode_func = build_ele_func(
            funcs, dfuncs,
            collect(funcs_vars_ntp[ele_input_names]),
            collect(funcs_vars_ntp[ele_state_names]),
            collect(funcs_params_ntp),
            funcs_nn_params
        )

        return new(
            name,
            nameinfo,
            funcs,
            dfuncs,
            flux_func,
            ode_func,
        )
    end
end

function build_ele_func(
    funcs::Vector{<:AbstractFlux},
    dfuncs::Vector{<:AbstractStateFlux},
    funcs_inputs::AbstractVector,
    funcs_states::AbstractVector,
    funcs_params::AbstractVector,
    funcs_nns::AbstractVector
)
    #* Call the method in SymbolicUtils.jl to build the Flux Function
    assign_list = Assignment[]
    #* get all output flux
    output_list = Num[]
    for func in funcs
        if func isa AbstractNeuralFlux
            #* For NeuralFlux, use nn_input to match the input variable
            push!(assign_list, Assignment(func.nn_info[:input], MakeArray(collect(func.input_info), Vector)))
            #* For NeuralFlux, use nn_output to match the calculation result of nn expr
            push!(assign_list, Assignment(func.nn_info[:output], func.flux_expr))
            #* According to the output variable name, match each index of nn_output
            for (idx, output) in enumerate(collect(func.output_info))
                push!(assign_list, Assignment(output, func.nn_info[:output][idx]))
                push!(output_list, output)
            end
        else
            #* According to the output variable name, match each result of the flux exprs
            for (output, expr) in zip(collect(func.output_info), func.flux_exprs)
                push!(assign_list, Assignment(output, expr))
                push!(output_list, output)
            end
        end
    end
    #* convert flux output to array
    flux_output_array = MakeArray(output_list, Vector)

    #* Set the input argument of ODE Function
    func_args = [
        #* argument 1: Function input and state variables
        DestructuredArgs(vcat(funcs_inputs, funcs_states)),
        #* argument 2: Function calculation parameters
        DestructuredArgs(funcs_params),
        #* argument 3: Function neuralnetwork parameters
        DestructuredArgs(funcs_nns),
    ]

    #* Construct Flux Function: Func(args, kwargs, body), where body represents the matching formula between each variable and expression
    merged_flux_func = @RuntimeGeneratedFunction(
        toexpr(Func(func_args, [], Let(assign_list, flux_output_array, false)))
    )

    if length(dfuncs) > 0
        #* convert diff state output to array
        diffst_output_array = MakeArray([dfunc.state_expr for dfunc in dfuncs], Vector)
        #* Set the input argument of ODE Function
        dfunc_args = [
            #* argument 1: Function input variables
            DestructuredArgs(funcs_inputs),
            #* argument 2: Function state variables
            DestructuredArgs(funcs_states),
            #* argument 3: Function calculation parameters
            DestructuredArgs(funcs_params),
            #* argument 4: Function neuralnetwork parameters
            DestructuredArgs(funcs_nns),
        ]
        merged_state_func = @RuntimeGeneratedFunction(
            toexpr(Func(dfunc_args, [], Let(assign_list, diffst_output_array, false)))
        )
    else
        merged_state_func = nothing
    end
    merged_flux_func, merged_state_func
end

function (ele::HydroElement)(
    input::AbstractMatrix,
    pas::ComponentVector;
    timeidx::Vector,
    solver::AbstractSolver=ODESolver(),
)
    # todo 输入的pas检验是否满足要求
    #* Extract the initial state of the parameters and element in the pas variable
    # todo When there is an ode function, it is necessary to solve the ode problem and obtain the state after the solution.
    if !isnothing(ele.ode_func)
        #* Call the solve_prob method to solve the state of element at the specified timeidx
        solved_states = solve_single_prob(ele, input=input, pas=pas, timeidx=timeidx, solver=solver)
        if solved_states == false
            solved_states = zeros(length(ele.nameinfo[:state]), length(timeidx))
        end
        #* Store the solved element state in fluxes
        fluxes = cat(input, solved_states, dims=1)
    else
        fluxes = input
        solved_states = nothing
    end
    #* excute other fluxes formula
    ele_output_matrix = run_fluxes(ele, input=fluxes, pas=pas)
    #* merge output and state
    if isnothing(solved_states)
        output_matrix = ele_output_matrix
    else
        output_matrix = cat(solved_states, ele_output_matrix, dims=1)
    end
    output_matrix
end

function run_fluxes(
    ele::HydroElement;
    #* var num * ts len
    input::AbstractMatrix,
    pas::ComponentVector,
)
    params_vec = collect([pas[:params][nm] for nm in ele.nameinfo[:param]])
    if length(ele.nameinfo[:nn]) > 0
        nn_params_vec = collect([pas[:nn][nm] for nm in ele.nameinfo[:nn]])
    else
        nn_params_vec = nothing
    end
    ele_output = ele.flux_func.(eachcol(input), Ref(params_vec), Ref(nn_params_vec))
    #* convert vector{vector} to matrix
    ele_output_matrix = reduce(hcat, ele_output)
    ele_output_matrix
end

function solve_single_prob(
    ele::HydroElement;
    input::AbstractMatrix,
    pas::ComponentVector,
    timeidx::Vector=collect(1:length(input[1])),
    solver::AbstractSolver=ODESolver(),
)
    params, init_states = pas[:params], pas[:initstates]

    #* Interpolate the input data. Since ordinary differential calculation is required, the data input must be continuous,
    #* so an interpolation function can be constructed to apply to each time point.
    itpfunc_list = map(eachrow(input)) do var
        LinearInterpolation(var, timeidx, extrapolate=true)
    end
    #* Construct a function for the ode_func input variable. Because of the difference in t, the ode_func input is not fixed.
    ode_input_func = (t) -> [itpfunc(t) for itpfunc in itpfunc_list]

    #* Extract the idx range of each variable in params,
    #* this extraction method is significantly more efficient than extracting by name
    params_idx = [getaxes(params)[1][nm].idx for nm in ele.nameinfo[:param]]
    # #* The input format is the input variable plus the intermediate state, which is consistent with the input of ode_func (see in line 45)
    ode_param_func = (p) -> [p[:params][idx] for idx in params_idx]

    if length(ele.nameinfo[:nn]) > 0
        nn_params = pas[:nn]
        nn_params_idx = [getaxes(nn_params)[1][nm].idx for nm in ele.nameinfo[:nn]]
        ode_nn_param_func = (p) -> [p[:nn][idx] for idx in nn_params_idx]
    else
        ode_nn_param_func = (_) -> nothing
    end

    #* Construct a temporary function that couples multiple ode functions to construct the solution for all states under the element
    function singel_ele_ode_func!(du, u, p, t)
        ode_input = ode_input_func(t)
        ode_params = ode_param_func(p)
        nn_params = ode_nn_param_func(p)
        du[:] = ele.ode_func(ode_input, u, ode_params, nn_params)
    end

    if typeof(solver) == ODESolver
        #* Construct ODEProblem based on DifferentialEquations.jl from this temporary function `singel_ele_ode_func!`
        prob = ODEProblem(
            singel_ele_ode_func!,
            collect(init_states[ele.nameinfo[:state]]),
            (timeidx[1], timeidx[end]),
            pas
        )
    elseif typeof(solver) == DiscreteSolver
        #* Construct ODEProblem based on DifferentialEquations.jl from this temporary function `singel_ele_ode_func!`
        prob = DiscreteProblem(
            singel_ele_ode_func!,
            collect(init_states[ele.nameinfo[:state]]),
            (timeidx[1], timeidx[end]),
            pas
        )
    end
    #* Solve the problem using the solver wrapper
    sol = solver(prob)
    sol
end

function run_multi_fluxes(
    ele::HydroElement;
    #* var num * node num * ts len
    input::AbstractArray,
    pas::ComponentVector,
    ptypes::Vector{Symbol}=collect(keys(pas[:params])),
)
    ele_param_vec = collect([collect([pas[:params][ptype][pname] for pname in ele.nameinfo[:param]]) for ptype in ptypes])
    if length(ele.nameinfo[:nn]) > 0
        nn_params_vec = collect([pas[:nn][nm] for nm in ele.nameinfo[:nn]])
    else
        nn_params_vec = nothing
    end
    #* array dims: (num of node, sequence length, variable dim)
    ele_output_vec = [ele.flux_func.(eachslice(input[:, i, :], dims=2), Ref(ele_param_vec[i]), Ref(nn_params_vec)) for i in 1:size(input)[2]]
    ele_output_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), [reduce(hcat, u) for u in ele_output_vec])
    ele_output_arr
end

function solve_multi_prob(
    ele::HydroElement;
    input::AbstractArray,
    pas::ComponentVector,
    timeidx::Vector,
    ptypes::Vector{Symbol}=collect(keys(pas[:params])),
    solver::LumpedHydro.AbstractSolver=LumpedHydro.ODESolver()
)
    #* 针对多个相同的state function采用并行化计算,这样能够避免神经网络反复多次计算减少梯度计算反馈
    #* 同时将多组state function放到同一个ode function中,这种并行计算或能提高预测性能,
    #* 这样每个步长的输入维度就是:节点个数*输入变量数
    #* 当前只针对unit相同的同步求解:
    params, init_states = pas[:params], pas[:initstates]
    # todo 保证两者长度一致
    # @assert length(ptypes) == length(input)
    
    #* Interpolate the input data. Since ordinary differential calculation is required, the data input must be continuous,
    #* so an interpolation function can be constructed to apply to each time point.
    itpfunc_vecs = [LinearInterpolation.(eachslice(i, dims=1), Ref(timeidx), extrapolate=true) for i in eachslice(input, dims=2)]
    ode_input_func = (t) -> [[itpfunc(t) for itpfunc in itpfunc_vec] for itpfunc_vec in itpfunc_vecs]
    
    #* 准备初始状态
    init_states_vec = collect([collect(init_states[ptype][ele.nameinfo[:state]]) for ptype in ptypes])
    init_states_matrix = reduce(hcat, init_states_vec)

    #* Extract the idx range of each variable in params,
    #* this extraction method is significantly more efficient than extracting by name
    params_idx = [getaxes(params[ptypes[1]])[1][nm].idx for nm in ele.nameinfo[:param]]
    #* Construct a function for the ode_func input variable. Because of the difference in t, the ode_func input is not fixed.
    #* The input format is the input variable plus the intermediate state, which is consistent with the input of ode_func (see in line 45)
    ode_param_func = (p) -> [[p[:params][ptype][idx] for idx in params_idx] for ptype in ptypes]

    #* 准备神经网络的参数
    if length(ele.nameinfo[:nn]) > 0
        nn_params_idx = [getaxes(pas[:nn])[1][nm].idx for nm in ele.nameinfo[:nn]]
        ode_nn_param_func = (p) -> [p[:nn][idx] for idx in nn_params_idx]
    else
        ode_nn_param_func = (_) -> nothing
    end

    #* Construct a temporary function that couples multiple ode functions to construct the solution for all states under the element
    function multi_ele_ode_func!(du, u, p, t)
        ode_input = ode_input_func(t)
        ode_params = ode_param_func(p)
        nn_params = ode_nn_param_func(p)
        tmp_output_vec = ele.ode_func.(ode_input, eachslice(u, dims=2), ode_params, Ref(nn_params))
        tmp_output = reduce(hcat, tmp_output_vec)
        du[:] = tmp_output
    end

    #* Construct ODEProblem based on DifferentialEquations.jl from this temporary function `single_ele_ode_func!`
    prob = ODEProblem(
        multi_ele_ode_func!,
        init_states_matrix,
        (timeidx[1], timeidx[end]),
        pas
    )

    #* Solve the problem using the solver wrapper
    solver(prob)
end


struct LagElement <: AbstractLagElement
    "the name of hydrological computation element "
    name::Symbol
    "element information: keys contains: input, output, param, state"
    nameinfo::NamedTuple
    "lag hydrological fluxes, used to flood routing"
    lfuncs::Vector

    function LagElement(
        name::Symbol;
        lfuncs::Vector{<:AbstractLagFlux},
    )
        #* Extract all variable names of funcs and dfuncs
        ele_input_names, ele_output_names = reduce(union, get_input_names.(lfuncs)), reduce(union, get_output_names.(lfuncs))
        #* Extract all parameters names of funcs and dfuncs
        ele_param_names = unique(reduce(union, get_param_names.(lfuncs)))
        #* Setup the name information of the hydroelement
        nameinfo = (input=ele_input_names, output=ele_output_names, param=ele_param_names)

        return new(
            name,
            nameinfo,
            lfuncs,
        )
    end
end

function (ele::LagElement)(
    input::AbstractMatrix,
    pas::ComponentVector;
    timeidx::Vector,
    solver::AbstractSolver
)
    #* Extract the initial state of the parameters and element in the pas variable
    params = pas[:params]
    lag_weights = [lfunc.lag_func(params[get_param_names(lfunc)[1]]) for lfunc in ele.lfuncs]

    function solve_lag_flux(input_vec, lag_weight)
        #* 首先将lagflux转换为discrete problem
        function lag_prob(u, p, t)
            u = circshift(u, -1)
            u[end] = 0.0
            input_vec[Int(t)] .* p[:weight] .+ u
        end

        prob = DiscreteProblem(lag_prob, lag_weight, (timeidx[1], timeidx[end]),
            ComponentVector(weight=lag_weight))
        #* 求解这个问题
        sol = solve(prob, FunctionMap())
        sol
    end

    sols = solve_lag_flux.(eachslice(input, dims=1), lag_weights)
    reduce(hcat, [Array(sol)[1, :] for sol in sols])'
end

function run_multi_fluxes(
    ele::LagElement;
    input::AbstractArray,
    pas::ComponentVector,
    ptypes::Vector{Symbol}=collect(keys(pas[:params])),
)
    #* array dims: (variable dim, num of node, sequence length)
    #* Extract the initial state of the parameters and element in the pas variable
    #* var_name * [weight_len * node_num]
    params = pas[:params]
    lag_weights = [[lfunc.lag_func(params[ptype][get_param_names(lfunc)[1]]) for ptype in ptypes] for lfunc in ele.lfuncs]

    function solve_lag_flux(input_vec, lag_weight)
        #* 首先将lagflux转换为discrete problem
        function lag_prob(u, p, t)
            tmp_u = circshift(u, -1)
            tmp_u[end] = 0.0
            input_vec[Int(t)] .* p[:w] .+ tmp_u
        end
        prob = DiscreteProblem(lag_prob, lag_weight, (1, length(input_vec)), ComponentVector(w=lag_weight))
        sol = solve(prob, FunctionMap())
        Array(sol)[1, :]
    end

    sols = map(eachindex(ele.lfuncs)) do idx
        node_sols = reduce(hcat, solve_lag_flux.(eachslice(input[idx, :, :], dims=1), lag_weights[idx]))
        node_sols
    end
    if length(sols) > 1
        sol_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), sols)
        return permutedims(sol_arr, (3, 1, 2))
    else
        return reshape(sols[1], 1, size(input)[3], size(input)[2])
    end
end
