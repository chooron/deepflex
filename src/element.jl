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
    The calculation topology map constructed based on common hydrological fluxes,
    ensures the orderly calculation of multiple hydrological fluxes.
    """
    topology::MetaTopology
    """
    Modelingtoolkit.jl related variables,
    This is a pre-built ODESystem based on the input hydrological flux,
    which is used to support the construction of the ODESystem after data input.
    """
    system::Union{ODESystem,Nothing}

    function HydroElement(
        name::Symbol;
        funcs::Vector,
        dfuncs::Vector=SimpleFlux[],
        mtk::Bool=true
    )
        topology = build_compute_topology(funcs)
        if !mtk
            return new{mtk}(
                name,
                funcs,
                dfuncs,
                topology,
                nothing
            )
        else
            var_names = unique(vcat(get_var_names(funcs, dfuncs)...))
            param_names = get_param_names(vcat(funcs, dfuncs))

            system = build_ele_system(funcs, dfuncs, var_names=var_names, param_names=param_names, name=name)

            return new{mtk}(
                name,
                funcs,
                dfuncs,
                topology,
                system
            )
        end
    end
end

# 求解并计算
function (ele::HydroElement)(
    input::NamedTuple,
    params::ComponentVector;
)
    fluxes = input
    for flux_idx in topological_sort(ele.topology.digraph)
        tmp_flux_name = ele.topology.node_names[flux_idx]
        if !(tmp_flux_name in keys(fluxes))
            tmp_func = ele.topology.node_maps[tmp_flux_name]
            tmp_output = tmp_func(fluxes, params)
            fluxes = merge(fluxes, tmp_output)
        end
    end
    fluxes
end

function (elements::AbstractVector{<:HydroElement})(
    input::NamedTuple,
    params::ComponentVector,
    topology::MetaTopology,
)
    #* 针对多个element的网络计算
    fluxes = input
    for flux_idx in topological_sort(topology.digraph)
        tmp_flux_name = topology.node_names[flux_idx]
        if !(tmp_flux_name in keys(fluxes))
            tmp_ele = topology.node_maps[tmp_flux_name]
            fluxes = merge(fluxes, tmp_ele(fluxes, params))
        end
    end
    fluxes
end

function add_inputflux!(
    ele::HydroElement;
    func_ntp::NamedTuple,
)
    for key in keys(func_ntp)
        for func in func_ntp[key]
            push!(ele.funcs, func)
        end
        dfunc = first(filter!(dfn -> dfn.output_names == key, ele.dfuncs))
        for func in func_ntp[key]
            dfunc.influx_names = vcat(dfunc.influx_names, get_output_names(func))
        end
    end
end

function add_outputflux!(
    ele::HydroElement;
    func_ntp::NamedTuple,
)
    for key in keys(func_ntp)
        for func in func_ntp[key]
            push!(ele.funcs, func)
        end
        dfunc = first(filter!(dfn -> dfn.output_names == key, ele.dfuncs))
        for func in func_ntp[key]
            dfunc.outflux_names = vcat(dfunc.outflux_names, get_output_names(func))
        end
    end
end

function solve_prob(
    ele::HydroElement{true};
    input::NamedTuple,
    params::ComponentVector,
    init_states::ComponentVector,
    solver::AbstractSolver=ODESolver()
)
    #* 准备计算数据，包括时间，名称，初始状态，neuralflux信息
    ts = collect(input[:time])
    ele_input_names, ele_state_names = get_input_names(ele), get_state_names(ele)
    nfunc_ntp = extract_neuralflux_ntp(ele.funcs)
    input_0 = namedtuple(ele_input_names, [input[nm][1] for nm in ele_input_names])
    initstates_ntp = namedtuple(keys(init_states), [init_states[nm] for nm in keys(init_states)])
    sol_0 = ele(merge(input_0, initstates_ntp), params)

    # 构建问题和setter
    build_sys = setup_input(ele.system, input=input[ele_input_names], time=ts, name=ele.name)
    prob = init_prob(build_sys, [ele.system], [nfunc_ntp], ts)

    # 获取u0和p实际值
    u0 = get_mtk_initstates([ele.system], sol_0, [ele_state_names], [nfunc_ntp])
    p = get_mtk_params([ele.system], params)

    new_prob = remake(prob, p=p, u0=u0)

    solved_states = solver(new_prob, ele_state_names)
    solved_states
end

function solve_prob(
    ele::HydroElement{false};
    input::NamedTuple,
    params::ComponentVector,
    init_states::ComponentVector,
    solver::AbstractSolver=ODESolver()
)
    ele_input_names = get_input_names(ele)
    ele_state_names = get_state_names(ele)
    itp_dict = namedtuple(ele_input_names, [LinearInterpolation(input[nm], input[:time], extrapolate=true) for nm in ele_input_names])

    function singel_ele_ode_func!(du, u, p, t)
        tmp_input = namedtuple(ele_input_names, [itp_dict[nm](t) for nm in ele_input_names])
        tmp_states = namedtuple(ele_state_names, u)
        tmp_fluxes = merge(tmp_input, tmp_states)
        tmp_fluxes = ele(tmp_fluxes, params)
        du[:] = [dfunc(tmp_fluxes, params)[dfunc.state_names] for dfunc in ele.dfuncs]
    end
    
    tspan = (input[:time][1], input[:time][end])
    prob = ODEProblem(singel_ele_ode_func!, collect(init_states[ele_state_names]), tspan)
    solved_states = solver(prob, ele_state_names)
    solved_states
end
