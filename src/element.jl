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
struct HydroElement <: AbstractElement
    "the name of hydrological computation element "
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
        #* Setup the name information of the hydroelement
        nameinfo = (input=ele_input_names, output=ele_output_names, state=ele_state_names, param=ele_param_names)
        #* Construct a NamedTuple of all func's variables and parameters
        funcs_vars_ntp = reduce(merge, [merge(func.input_info, func.output_info) for func in vcat(funcs, dfuncs)])
        funcs_params_ntp = reduce(merge, [func.param_info for func in vcat(funcs, dfuncs)])
        #* Construct a function for ordinary differential calculation based on dfunc and funcs
        flux_func, ode_func = build_ele_func(funcs, dfuncs, funcs_vars_ntp[ele_input_names], funcs_vars_ntp[ele_state_names], funcs_params_ntp)

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
    funcs_input_ntp::NamedTuple,
    funcs_state_ntp::NamedTuple,
    funcs_params_ntp::NamedTuple
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
        DestructuredArgs(vcat(collect(funcs_input_ntp), collect(funcs_state_ntp))),
        #* argument 2: Function calculation parameters
        DestructuredArgs(collect(funcs_params_ntp))
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
            DestructuredArgs(collect(funcs_input_ntp)),
            #* argument 2: Function state variables
            DestructuredArgs(collect(funcs_state_ntp)),
            #* argument 3: Function calculation parameters
            DestructuredArgs(collect(funcs_params_ntp))
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
    #* Extract the initial state of the parameters and element in the pas variable
    params, init_states = pas[:params], pas[:initstates]
    # todo When there is an ode function, it is necessary to solve the ode problem and obtain the state after the solution.
    if !isnothing(ele.ode_func)
        #* Call the solve_prob method to solve the state of element at the specified timeidx
        solved_states = solve_single_prob(ele, input=input, params=params, init_states=init_states, timeidx=timeidx, solver=solver)
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
    ele_params_vec = collect([params[nm] for nm in ele.nameinfo[:param]])
    ele_output = ele.flux_func.(eachcol(fluxes), Ref(ele_params_vec))
    #* convert vector{vector} to matrix
    ele_output_matrix = reduce(hcat, ele_output)
    #* merge output and state
    if isnothing(solved_states)
        output_matrix = ele_output_matrix
    else
        output_matrix = cat(solved_states, ele_output_matrix, dims=1)
    end
    output_matrix
end

function solve_single_prob(
    ele::HydroElement;
    input::AbstractMatrix,
    params::ComponentVector,
    init_states::ComponentVector,
    timeidx::Vector=collect(1:length(input[1])),
    solver::AbstractSolver=ODESolver(),
)
    #* Extract the idx range of each variable in params,
    #* this extraction method is significantly more efficient than extracting by name
    params_idx = [getaxes(params)[1][nm].idx for nm in ele.nameinfo[:param]]

    #* Interpolate the input data. Since ordinary differential calculation is required, the data input must be continuous,
    #* so an interpolation function can be constructed to apply to each time point.
    itpfunc_list = map(eachrow(input)) do var
        LinearInterpolation(var, timeidx, extrapolate=true)
    end

    #* Construct a function for the ode_func input variable. Because of the difference in t, the ode_func input is not fixed.
    #* The input format is the input variable plus the intermediate state, which is consistent with the input of ode_func (see in line 45)
    ode_input_func = (t) -> [itpfunc(t) for itpfunc in itpfunc_list]
    ode_param_func = (p) -> [p[idx] for idx in params_idx]

    if typeof(solver) == ODESolver
        #* Construct a temporary function that couples multiple ode functions to construct the solution for all states under the element
        function singel_ele_ode_func!(du, u, p, t)
            du[:] = ele.ode_func(ode_input_func(t), u, ode_param_func(p))
        end

        #* Construct ODEProblem based on DifferentialEquations.jl from this temporary function `singel_ele_ode_func!`
        prob = ODEProblem(
            singel_ele_ode_func!,
            collect(init_states[ele.nameinfo[:state]]),
            (timeidx[1], timeidx[end]),
            params
        )
    elseif typeof(solver) == DiscreteSolver
        function singel_ele_disc_func!(du, u, p, t)
            du[:] = ele.ode_func(ode_input_func(t, u), ode_param_func(p))
        end

        #* Construct ODEProblem based on DifferentialEquations.jl from this temporary function `singel_ele_ode_func!`
        prob = DiscreteProblem(
            singel_ele_disc_func!,
            collect(init_states[ele_state_names]),
            (timeidx[1], timeidx[end]),
            params
        )
    end
    #* Solve the problem using the solver wrapper
    sol = solver(prob)
    sol
end

function run_multi_fluxes(
    ele::HydroElement;
    #* node num * ts len * var num
    input::AbstractArray,
    params::ComponentVector,
)
    ele_param_vec = collect([collect([params[n][nn] for nn in ele.nameinfo[:param]]) for n in keys(params)])
    #* array dims: (num of node, sequence length, variable dim)
    ele_output_vec = [ele.flux_func.(eachslice(input[:, i, :], dims=2), Ref(ele_param_vec[i])) for i in 1:size(input)[2]]
    ele_output_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), [reduce(hcat, u) for u in ele_output_vec])
    ele_output_arr
end


function solve_multi_prob(
    ele::HydroElement;
    input::AbstractArray,
    params::ComponentVector,
    init_states::ComponentVector,
    timeidx::Vector,
    solver::LumpedHydro.AbstractSolver=LumpedHydro.ODESolver()
)
    #* 针对多个相同的state function采用并行化计算,这样能够避免神经网络反复多次计算减少梯度计算反馈
    #* 同时将多组state function放到同一个ode function中,这种并行计算或能提高预测性能,
    #* 这样每个步长的输入维度就是:节点个数*输入变量数
    #* 当前只针对unit相同的同步求解:
    #* Extract the idx range of each variable in params,
    #* this extraction method is significantly more efficient than extracting by name
    param_keys = keys(params)
    params_idx = [getaxes(params[param_keys[1]])[1][nm].idx for nm in ele.nameinfo[:param]]

    #* 准备初始状态
    init_states_vec = collect([collect(init_states[ndname][ele.nameinfo[:state]]) for ndname in keys(init_states)])
    init_states_matrix = reduce(hcat, init_states_vec)

    #* Interpolate the input data. Since ordinary differential calculation is required, the data input must be continuous,
    #* so an interpolation function can be constructed to apply to each time point.
    itpfunc_vecs = [LinearInterpolation.(eachslice(i, dims=1), Ref(timeidx), extrapolate=true) for i in eachslice(input, dims=2)]

    #* Construct a function for the ode_func input variable. Because of the difference in t, the ode_func input is not fixed.
    #* The input format is the input variable plus the intermediate state, which is consistent with the input of ode_func (see in line 45)
    ode_input_func = (t) -> [[itpfunc(t) for itpfunc in itpfunc_vec] for itpfunc_vec in itpfunc_vecs]
    ode_param_func = (p) -> [[p[pname][idx] for idx in params_idx] for pname in param_keys]

    #* Construct a temporary function that couples multiple ode functions to construct the solution for all states under the element
    function multi_ele_ode_func!(du, u, p, t)
        ode_input = ode_input_func(t)
        ode_params = ode_param_func(p)
        tmp_output_vec = ele.ode_func.(ode_input, eachslice(u, dims=2), ode_params)
        tmp_output = reduce(hcat, tmp_output_vec)
        du[:] = tmp_output
    end

    #* Construct ODEProblem based on DifferentialEquations.jl from this temporary function `single_ele_ode_func!`
    prob = ODEProblem(
        multi_ele_ode_func!,
        init_states_matrix,
        (timeidx[1], timeidx[end]),
        params
    )

    #* Solve the problem using the solver wrapper
    solver(prob)
end
