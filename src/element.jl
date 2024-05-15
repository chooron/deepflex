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
    topology::SimpleDiGraph
    """
    Modelingtoolkit.jl related variables,
    define the variables according to the input and output names of the hydrological flux,
    and save them as NamedTuple
    """
    varinfo::NamedTuple
    """
    Modelingtoolkit.jl related variables,
    define the parameters according to the parameter names of the hydrological flux,
    and save them as NamedTuple
    """
    paraminfo::NamedTuple
    """
    Modelingtoolkit.jl related variables,
    This is a pre-built ODESystem based on the input hydrological flux,
    which is used to support the construction of the ODESystem after data input.
    """
    sys::Union{ODESystem,Nothing}

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
                NamedTuple(),
                NamedTuple(),
                nothing
            )
        else
            input_names = get_input_names(funcs, dfuncs)
            output_names = get_output_names(funcs)
            state_names = get_state_names(dfuncs)
            param_names = get_param_names(vcat(funcs, dfuncs))

            varinfo, paraminfo = init_var_param(input_names, output_names, state_names, param_names)
            sys = build_ele_system(funcs, dfuncs, varinfo, paraminfo, name=name)

            return new{mtk}(
                name,
                funcs,
                dfuncs,
                topology,
                varinfo,
                paraminfo,
                sys
            )
        end
    end
end

# 求解并计算
function (ele::HydroElement)(
    input::NamedTuple,
    pas::ComponentVector;
)
    #* 构建各个输出变量对应的计算函数
    func_ntp = namedtuple(
        vcat([get_output_names(func) for func in ele.funcs]...),
        vcat([repeat([func], length(get_output_names(func))) for func in ele.funcs]...)
    )
    ele_var_names = unique(vcat(get_input_output_names(ele.funcs)...))
    fluxes = input
    for flux_idx in topological_sort(ele.topology)
        tmp_flux_name = ele_var_names[flux_idx]
        if !(tmp_flux_name in keys(fluxes))
            tmp_func = func_ntp[tmp_flux_name]
            tmp_output = tmp_func(fluxes, pas[:params])
            fluxes = merge(fluxes, tmp_output)
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

function get_sol_0(
    ele::HydroElement{true};
    input::NamedTuple,
    params::ComponentVector,
    init_states::ComponentVector
)
    ele_input_names = setdiff(get_input_names(ele.funcs), keys(init_states))
    input0 = namedtuple(ele_input_names, [input[nm][1] for nm in ele_input_names])
    init_states_ntp = namedtuple(keys(init_states), [init_states[nm] for nm in keys(init_states)])
    sol_0 = merge(input0, init_states_ntp)
    for func in ele.funcs
        sol_0 = merge(sol_0, func(sol_0, params))
    end
    ComponentVector(sol_0)
end

function setup_prob(
    ele::HydroElement{true},
    prob::ODEProblem;
    input::NamedTuple,
    params::ComponentVector,
    init_states::ComponentVector,
    nfunc_ntp::NamedTuple
)
    sol_0 = get_sol_0(ele, input=input, params=params, init_states=init_states)
    u0 = get_mtk_initstates(ele.sys; sol_0=sol_0, state_names=get_state_names(ele), nfunc_ntp=nfunc_ntp)
    p = get_mtk_params(ele.sys; params=params)
    remake(prob, p=p, u0=u0)
end

function solve_prob(
    ele::HydroElement{true};
    input::NamedTuple,
    pas::ComponentVector,
    solver::AbstractSolver=ODESolver()
)
    params = pas[:params] isa Vector ? ComponentVector() : pas[:params]
    init_states = pas[:initstates] isa Vector ? ComponentVector() : pas[:initstates]

    ele_input_names = get_input_names(ele)
    nfunc_ntp = extract_neuralflux_ntp(ele.funcs)

    prob = setup_input(ele.sys, input=input[ele_input_names], time=input[:time], nfunc_ntp=nfunc_ntp, name=ele.name)
    new_prob = setup_prob(ele, prob, input=input, params=params, init_states=init_states, nfunc_ntp=nfunc_ntp)

    solved_states = solver(new_prob, get_state_names(ele.dfuncs))
    solved_states
end

function solve_prob(
    ele::HydroElement{false};
    input::NamedTuple,
    pas::Union{NamedTuple,ComponentVector},
    solver::AbstractSolver=ODESolver()
)
    params = pas[:params] isa Vector ? ComponentVector() : pas[:params]
    init_states = pas[:initstates] isa Vector ? ComponentVector() : pas[:initstates]
    ele_input_names = get_input_names(ele)
    itp_dict = namedtuple(ele_input_names, [LinearInterpolation(input[nm], input[:time], extrapolate=true) for nm in ele_input_names])
    function singel_ele_ode_func!(du, u, p, t)
        tmp_input = namedtuple(ele_input_names, [itp_dict[nm](t) for nm in ele_input_names])
        tmp_states = namedtuple(get_state_names(ele), u)
        tmp_fluxes = merge(tmp_input, tmp_states)
        for func in ele.funcs
            tmp_fluxes = merge(tmp_fluxes, func(tmp_fluxes, params))
        end
        du[:] = [dfunc(tmp_fluxes, params)[first(get_output_names(dfunc))] for dfunc in ele.dfuncs]
    end
    prob = ODEProblem(singel_ele_ode_func!, collect(init_states[get_state_names(ele)]), (input[:time][1], input[:time][end]))
    solved_states = solver(prob, get_state_names(ele))
    solved_states
end

export HydroElement, add_inputflux!, add_outputflux!, solve_prob