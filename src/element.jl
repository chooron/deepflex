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

    function HydroElement(
        name::Symbol;
        funcs::Vector,
        dfuncs::Vector=StateFlux[],
        lfuncs::Vector=LagFlux[],
    )
        ele_input_names, ele_output_names, ele_state_names = get_var_names(funcs, dfuncs)
        ele_param_names = get_param_names(vcat(funcs, dfuncs))
        nameinfo = (input=ele_input_names, output=ele_output_names, state=ele_state_names, param=ele_param_names)

        return new(
            name,
            nameinfo,
            funcs,
            dfuncs,
            lfuncs,
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
    params = pas[:params]
    init_states = pas[:initstates]
    ele_state_names = ele.nameinfo[:state]
    if !solved & (length(ele.dfuncs) > 0)
        solved_states = solve_prob(ele, input=fluxes, params=params, init_states=init_states, timeidx=timeidx, solver=solver)
        solved_states = ComponentVector(NamedTuple{Tuple(ele_state_names)}(solved_states))
        fluxes = merge(fluxes, solved_states)
    end
    for ele_func in ele.funcs
        tmp_mtr = reduce(hcat, collect(fluxes[get_input_names(ele_func)]))
        tmp_params_vec = eltype(params)[params[get_param_names(ele_func)]...]
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
    params::ComponentVector,
    init_states::ComponentVector,
    timeidx::Vector=collect(1:length(input[1])),
    solver::AbstractSolver=ODESolver(saveat=timeidx),
)
    ele_input_names = ele.nameinfo[:input]
    ele_state_names = ele.nameinfo[:state]

    itpfunc_ntp = NamedTuple{Tuple(ele_input_names)}(
        [LinearInterpolation(input[nm], timeidx, extrapolate=true) for nm in ele_input_names]
    )

    dfuncs_input_names_list = [get_input_names(dfunc) for dfunc in ele.dfuncs]

    dfuncs_input_func_list = [
        (t, u) -> begin
            [
                begin
                    state_idx = findfirst(x -> x == nm, ele_state_names)
                    isnothing(state_idx) ? itpfunc_ntp[nm](t) : u[state_idx]
                end
                for nm in dfuncs_input_names
            ]
        end
        for dfuncs_input_names in dfuncs_input_names_list
    ]

    if solver isa ODESolver
        function singel_ele_ode_func!(du, u, p, t)
            du[:] = [dfunc.inner_func(dfuncs_input_func(t, u), p) for (dfunc, dfuncs_input_func) in zip(ele.dfuncs, dfuncs_input_func_list)]
        end

        prob = ODEProblem(
            singel_ele_ode_func!,
            collect(init_states[ele_state_names]),
            (timeidx[1], timeidx[end]),
            [params[nm] for nm in ele.nameinfo[:param]]
        )
    elseif solver isa DiscreteSolver
        function singel_ele_ode_func(u, p, t)
            u[:] = [u[idx] + ele.dfuncs[idx].inner_func(dfuncs_input_func_list[idx](t, u), p)[idx] for idx in 1:length(ele.dfuncs)]
        end

        prob = DiscreteProblem(
            singel_ele_ode_func,
            collect(init_states[ele_state_names]),
            (timeidx[1], timeidx[end]),
            [params[nm] for nm in ele.nameinfo[:param]]
        )
    end

    solve_u = solver(prob)
    solve_u
end