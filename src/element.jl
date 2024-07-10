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
    "hydrological computation element name"
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
    Hydrological lag fluxes
    """
    ode_funcs::Vector

    function HydroElement(
        name::Symbol;
        funcs::Vector,
        dfuncs::Vector=StateFlux[],
        lfuncs::Vector=LagFlux[],
    )
        ele_input_names, ele_output_names, ele_state_names = get_var_names(funcs, dfuncs)
        ele_param_names = get_param_names(vcat(funcs, dfuncs))
        nameinfo = (input=ele_input_names, output=ele_output_names, state=ele_state_names, param=ele_param_names)
        ode_funcs = [build_state_func(funcs, dfunc, vcat(ele_input_names, ele_state_names)) for dfunc in dfuncs]

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
    fluxes::NamedTuple,
    pas::ComponentVector;
    timeidx::Vector=collect(1:length(fluxes[first(keys(fluxes))])),
    solver::AbstractSolver=ODESolver(saveat=timeidx),
    solved::Bool=false
)
    params = NamedTuple(pas[:params])
    init_states = NamedTuple(pas[:initstates])
    ele_state_names = ele.nameinfo[:state]
    if !solved & (length(ele.dfuncs) > 0)
        solved_states = solve_prob(ele, input=fluxes, params=params, init_states=init_states, timeidx=timeidx, solver=solver)
        solved_states = ComponentVector(NamedTuple{Tuple(ele_state_names)}(solved_states))
        fluxes = merge(fluxes, solved_states)
    end
    for ele_func in ele.funcs
        tmp_mtr = reduce(hcat, collect(fluxes[get_input_names(ele_func)]))
        tmp_params_vec = collect([params[nm] for nm in get_param_names(ele_func)])
        tmp_output = ele_func(tmp_mtr, tmp_params_vec)
        tmp_output_ntp = NamedTuple{Tuple(get_output_names(ele_func))}(eachcol(tmp_output))
        fluxes = merge(fluxes, tmp_output_ntp)
    end
    for ele_func in ele.lfuncs
        tmp_var_vec = fluxes[get_input_names(ele_func)[1]]
        tmp_lag_time = eltype(params)(params[get_param_names(ele_func)[1]])
        tmp_output = ele_func(tmp_var_vec, [tmp_lag_time])
        tmp_output_ntp = NamedTuple{Tuple(get_output_names(ele_func))}(tmp_output)
        fluxes = merge(fluxes, tmp_output_ntp)
    end
    fluxes
end

function solve_prob(
    ele::HydroElement;
    input::NamedTuple,
    params::NamedTuple,
    init_states::NamedTuple,
    timeidx::Vector=collect(1:length(input[1])),
    solver::AbstractSolver=ODESolver(saveat=timeidx),
)
    ele_input_names = ele.nameinfo[:input]
    ele_state_names = ele.nameinfo[:state]

    params = ComponentVector(params)
    params_idx = [getaxes(params)[1][nm].idx for nm in ele.nameinfo[:param]]

    ptypes = [Vector{eltype(params)}, Vector{eltype(params)}]

    itpfunc_ntp = NamedTuple{Tuple(ele_input_names)}(
        [LinearInterpolation(input[nm], timeidx, extrapolate=true) for nm in ele_input_names]
    )

    ode_input_func = (t, u) -> vcat([itpfunc_ntp[nm](t) for nm in ele_input_names], u)

    function singel_ele_ode_func!(du, u, p, t)
        du[:] = [ode_func(ode_input_func(t, u), [p[idx] for idx in params_idx], ptypes) for ode_func in ele.ode_funcs]
    end

    prob = ODEProblem(
        singel_ele_ode_func!,
        collect(init_states[ele_state_names]),
        (timeidx[1], timeidx[end]),
        params
    )

    solve_u = solver(prob)
    solve_u
end


function build_state_func(
    funcs::Vector{<:AbstractFlux},
    dfunc::AbstractStateFlux,
    input_names::Vector{Symbol},
)
    fluxes_vars_ntp = reduce(merge, [merge(func.input_info, func.output_info) for func in vcat(funcs, [dfunc])])
    funcs_params_ntp = reduce(merge, [func.param_info for func in vcat(funcs, [dfunc])])

    assign_list = Assignment[]
    for func in funcs
        if func isa AbstractNeuralFlux
            push!(assign_list, Assignment(func.nn_info[:input], MakeArray(collect(func.input_info), Vector)))
            push!(assign_list, Assignment(func.nn_info[:output], func.flux_expr))
            for (idx, output) in enumerate(collect(func.output_info))
                push!(assign_list, Assignment(output, func.nn_info[:output][idx]))
            end
        else
            for (output, expr) in zip(collect(func.output_info), func.flux_exprs)
                push!(assign_list, Assignment(output, expr))
            end
        end
    end

    func_args = [
        DestructuredArgs(collect(fluxes_vars_ntp[input_names])),
        DestructuredArgs(collect(funcs_params_ntp)),
        DestructuredArgs([func.nn_info[:chain_ptype] for func in funcs if func isa AbstractNeuralFlux])
    ]

    merged_state_func = @RuntimeGeneratedFunction(
        toexpr(Func(func_args, [], Let(assign_list, dfunc.state_expr, false)))
    )
    merged_state_func
end
