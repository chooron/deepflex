"""
$(TYPEDEF)
The basic hydrological calculation module contains multiple hydrological fluxes,
and can simulate the balance calculation of a physical module.
# Fields
$(FIELDS)
# Example
```
funcs = [
    PetFlux([:temp, :lday]),
    SnowfallFlux([:prcp, :temp], param_names=[:Tmin]),
    MeltFlux([:snowwater, :temp], param_names=[:Tmax, :Df]),
    RainfallFlux([:prcp, :temp], param_names=[:Tmin]),
    InfiltrationFlux([:rainfall, :melt])
]

dfuncs = [
    StateFlux([:snowfall], [:melt], :snowwater),
]

HydroElement(
    Symbol(name, :_surface_),
    funcs=funcs,
    dfuncs=dfuncs,
    mtk=mtk,
)
```
"""
struct HydroElement{mtk} <: AbstractElement
    "hydrological computation element name"
    name::Symbol
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
    Modelingtoolkit.jl related variables,
    This is a pre-built ODESystem based on the input hydrological flux,
    which is used to support the construction of the ODESystem after data input.
    """
    system::Union{DiscreteSystem,ODESystem,Nothing}

    function HydroElement(
        name::Symbol;
        funcs::Vector,
        dfuncs::Vector=StateFlux[],
        lfuncs::Vector=LagFlux[],
        mtk::Bool=false
    )
        return new{mtk}(
            name,
            funcs,
            dfuncs,
            lfuncs,
            mtk ? build_ele_system(funcs, dfuncs, name=name) : nothing
        )
    end
end

function (ele::HydroElement)(
    input::NamedTuple,
    pas::Union{ComponentVector,NamedTuple};
    timeidx::Vector=collect(1:length(input)),
    solver::AbstractSolver=ODESolver(),
    solved::Bool=false
)
    params = NamedTuple(pas[:params])
    init_states = NamedTuple(pas[:initstates])
    fluxes = input
    if !solved & (length(ele.dfuncs) > 0)
        solved_states = solve_prob(ele, input=input, params=params, init_states=init_states, timeidx=timeidx, solver=solver)
        solved_states = NamedTuple{Tuple(get_state_names(ele))}(solved_states)
        fluxes = merge(fluxes, solved_states)
    end
    for ele_func in vcat(ele.funcs, ele.lfuncs)
        tmp_mtr = hcat([fluxes[nm] for nm in get_input_names(ele_func)]...)
        tmp_params_vec = eltype(params)[params[get_param_names(ele_func)]...]
        tmp_output = ele_func(tmp_mtr, tmp_params_vec)
        tmp_output_ntp = NamedTuple{Tuple(get_output_names(ele_func))}(tmp_output)
        fluxes = merge(fluxes, tmp_output_ntp)
    end
    fluxes
end

function solve_prob(
    ele::HydroElement{true};
    input::NamedTuple,
    params::NamedTuple,
    init_states::NamedTuple,
    timeidx::Vector=collect(1:length(input[1])),
    solver::AbstractSolver=ODESolver()
)
    #* 准备计算数据，包括时间，名称，初始状态，neuralflux信息
    ele_input_names, _, ele_state_names = get_var_names(ele)
    ele_input_state_names = vcat(ele_input_names, ele_state_names)

    nn_flux_list = filter(flux -> flux isa AbstractNeuralFlux, ele.funcs)
    flux_chain_names = [flux.chain_name for flux in nn_flux_list]
    nfunc_ntp = length(nn_flux_list) == 0 ? NamedTuple() : NamedTuple{Tuple(flux_chain_names)}(
        [get_input_names(flux) for flux in nn_flux_list]
    )

    #* 构建问题
    build_sys = setup_input(ele.system, input, timeidx, ele_input_names, ele.name)
    prob = init_prob(build_sys, ele.system, nfunc_ntp, timeidx)

    #* 计算element在初始时间的水文通量
    #* 当神经水文通量数量大于0时需要计算所有初始状态的水文通量
    if length(nfunc_ntp) > 0
        sol_0 = NamedTuple{Tuple(ele_input_state_names)}(
            [nm in ele_input_names ? input[nm][1] : init_states[nm] for nm in ele_input_state_names]
        )
        for tmp_func in ele.funcs
            tmp_output = tmp_func(sol_0, params)
            for (idx, nm) in enumerate(get_output_names(tmp_func))
                sol_0 = merge(sol_0, NamedTuple{tuple(nm)}(tmp_output[idx]))
            end
        end
    else
        sol_0 = init_states
    end

    #* 获取u0和p实际值
    u0 = get_mtk_initstates([ele.system], sol_0, [ele_state_names], [nfunc_ntp])
    p = get_mtk_params([ele.system], params)

    new_prob = remake(prob, p=p, u0=u0)

    solved_states = solver(new_prob)
    solved_states
end

function solve_prob(
    ele::HydroElement{false};
    input::NamedTuple,
    params::NamedTuple,
    init_states::NamedTuple,
    timeidx::Vector=collect(1:length(input[1])),
    solver::AbstractSolver=ODESolver(),
)
    ele_input_names = get_input_names(ele)
    ele_state_names = get_state_names(ele)

    funcs_input_names, funcs_output_names = get_input_output_names(ele.funcs)
    dfuncs_input_names = [union(funcs_input_names, setdiff(get_input_names(dfunc), funcs_output_names)) for dfunc in ele.dfuncs]

    itpfunc_ntp = NamedTuple{Tuple(ele_input_names)}(
        [LinearInterpolation(input[nm], timeidx, extrapolate=true) for nm in ele_input_names]
    )

    dfuncs_input_func_list = [
        (t, u) -> begin
            [
                begin
                    state_idx = findfirst(x -> x == nm, ele_state_names)
                    isnothing(state_idx) ? itpfunc_ntp[nm](t) : u[state_idx]
                end
                for nm in dfuncs_input_name
            ]
        end
        for dfuncs_input_name in dfuncs_input_names
    ]

    if solver isa ODESolver
        function singel_ele_ode_func!(du, u, p, t)
            du[:] = [dfunc(dfuncs_input_func(t, u), p) for (dfunc, dfuncs_input_func) in zip(ele.dfuncs, dfuncs_input_func_list)]
        end

        prob = ODEProblem(
            singel_ele_ode_func!,
            collect(init_states[ele_state_names]),
            (timeidx[1], timeidx[end]),
            collect(params[get_param_names(ele)])
        )
    elseif solver isa DiscreteSolver
        function singel_ele_ode_func(u, p, t)
            u[:] = [u[idx] + ele.dfuncs[idx](dfuncs_input_func_list[idx](t, u), p)[idx] for idx in 1:length(ele.dfuncs)]
        end

        prob = DiscreteProblem(
            singel_ele_ode_func,
            collect(init_states[ele_state_names]),
            (timeidx[1], timeidx[end]),
            collect(params[get_param_names(ele)])
        )
    end

    solve_u = solver(prob)
    solve_u
end