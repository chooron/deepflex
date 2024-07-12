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
    Hydrological lag fluxes
    """
    lfuncs::Vector
    """
    Hydrological ode functions
    """
    ode_funcs::Vector

    function HydroElement(
        name::Symbol;
        funcs::Vector,
        dfuncs::Vector=StateFlux[],
        lfuncs::Vector=LagFlux[],
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
        ode_funcs = [build_state_func(funcs, dfunc, funcs_vars_ntp[vcat(ele_input_names, ele_state_names)], funcs_params_ntp) for dfunc in dfuncs]

        return new(
            name,
            nameinfo,
            funcs,
            dfuncs,
            lfuncs,
            ode_funcs
        )
    end
end

function (ele::HydroElement)(
    input::NamedTuple,
    pas::ComponentVector;
    timeidx::Vector=collect(1:length(input[first(keys(input))])),
    solver::AbstractSolver=ODESolver(saveat=timeidx),
    solved::Bool=false
)
    #* Define the NamedTuple class to store the calculation results of each func and use them for subsequent calculations
    fluxes = input
    #* Extract the initial state of the parameters and element in the pas variable
    params, init_states = pas[:params], pas[:initstates]
    #* When there is an ode function, it is necessary to solve the ode problem and obtain the state after the solution.
    if !solved & (length(ele.ode_funcs) > 0)
        #* Call the solve_prob method to solve the state of element at the specified timeidx
        solved_states = solve_prob(ele, input=fluxes, params=params, init_states=init_states, timeidx=timeidx, solver=solver)
        #* Store the solved element state in fluxes
        fluxes = merge(fluxes, NamedTuple{Tuple(ele.nameinfo[:state])}(solved_states))
    end
    #* After solving the ode problem, calculate other funcs
    for ele_func in vcat(ele.funcs, ele.lfuncs)
        #* AbstractLagFlux and AbstractSimpleFlux have different solutions.
        if ele_func isa AbstractLagFlux
            #* Extracting variables for routing (Vector{sequence length})
            tmp_var_vec = fluxes[get_input_names(ele_func)[1]]
            #* Extracting parameter for routing
            tmp_lag_time = params[get_param_names(ele_func)[1]]
            #* Get routing result based on the unit hydrograph
            tmp_output = ele_func(tmp_var_vec, [tmp_lag_time])
            #* Create a NamedTuple with (routing variable and routing results)
            tmp_output_ntp = NamedTuple{Tuple(get_output_names(ele_func))}(tmp_output)
        else
            #* Extracting variables for calculating (Matrix{sequence length, variable dim})
            tmp_mtr = reduce(hcat, collect(fluxes[get_input_names(ele_func)]))
            #* Extracting parameter for calculating
            tmp_params_vec = collect([params[nm] for nm in get_param_names(ele_func)])
            #* Get calculating results
            tmp_output = ele_func(tmp_mtr, tmp_params_vec)
            #* Create a NamedTuple with (output variable and output results)
            tmp_output_ntp = NamedTuple{Tuple(get_output_names(ele_func))}(eachcol(tmp_output))
        end
        #* Save the calculation results and use them for subsequent calculations
        fluxes = merge(fluxes, tmp_output_ntp)
    end
    fluxes
end

function solve_prob(
    ele::HydroElement;
    input::NamedTuple,
    params::ComponentVector,
    init_states::ComponentVector,
    timeidx::Vector=collect(1:length(input[1])),
    solver::AbstractSolver=ODESolver(saveat=timeidx),
)
    ele_input_names = ele.nameinfo[:input]
    ele_state_names = ele.nameinfo[:state]

    #* Extract the idx range of each variable in params,
    #* this extraction method is significantly more efficient than extracting by name
    params_idx = [getaxes(params)[1][nm].idx for nm in ele.nameinfo[:param]]

    #* Interpolate the input data. Since ordinary differential calculation is required, the data input must be continuous,
    #* so an interpolation function can be constructed to apply to each time point.
    itpfunc_ntp = NamedTuple{Tuple(ele_input_names)}(
        [LinearInterpolation(input[nm], timeidx, extrapolate=true) for nm in ele_input_names]
    )

    #* Construct a function for the ode_func input variable. Because of the difference in t, the ode_func input is not fixed.
    #* The input format is the input variable plus the intermediate state, which is consistent with the input of ode_func (see in line 45)
    ode_input_func = (t, u) -> vcat([itpfunc_ntp[nm](t) for nm in ele_input_names], u)

    #* Construct a temporary function that couples multiple ode functions to construct the solution for all states under the element
    function singel_ele_ode_func!(du, u, p, t)
        du[:] = [ode_func(ode_input_func(t, u), [p[idx] for idx in params_idx]) for ode_func in ele.ode_funcs]
    end

    #* Construct ODEProblem based on DifferentialEquations.jl from this temporary function `singel_ele_ode_func!`
    prob = ODEProblem(
        singel_ele_ode_func!,
        collect(init_states[ele_state_names]),
        (timeidx[1], timeidx[end]),
        params
    )

    #* Solve the problem using the solver wrapper
    solve_u = solver(prob)
    solve_u
end

function build_state_func(
    funcs::Vector{<:AbstractFlux},
    dfunc::AbstractStateFlux,
    funcs_input_ntp::NamedTuple,
    funcs_params_ntp::NamedTuple
)
    #* Call the method in SymbolicUtils.jl to build the ODE Function
    assign_list = Assignment[]
    for func in funcs
        if func isa AbstractNeuralFlux
            #* For NeuralFlux, use nn_input to match the input variable
            push!(assign_list, Assignment(func.nn_info[:input], MakeArray(collect(func.input_info), Vector)))
            #* For NeuralFlux, use nn_output to match the calculation result of nn expr
            push!(assign_list, Assignment(func.nn_info[:output], func.flux_expr))
            #* According to the output variable name, match each index of nn_output
            for (idx, output) in enumerate(collect(func.output_info))
                push!(assign_list, Assignment(output, func.nn_info[:output][idx]))
            end
        else
            #* According to the output variable name, match each result of the flux exprs
            for (output, expr) in zip(collect(func.output_info), func.flux_exprs)
                push!(assign_list, Assignment(output, expr))
            end
        end
    end
    #* Set the input argument of ODE Function
    func_args = [
        #* argument 1: Function input variables
        DestructuredArgs(collect(funcs_input_ntp)),
        #* argument 2: Function calculation parameters
        DestructuredArgs(collect(funcs_params_ntp))
    ]

    #* Construct ODE Function: Func(args, kwargs, body), where body represents the matching formula between each variable and expression
    merged_state_func = @RuntimeGeneratedFunction(
        toexpr(Func(func_args, [], Let(assign_list, dfunc.state_expr, false)))
    )
    merged_state_func
end
