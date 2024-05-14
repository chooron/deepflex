"""
HydroElement
"""
struct HydroElement{mtk} <: AbstractElement
    #* component名称
    name::Symbol
    #* fluxes and dfluxes
    funcs::Vector
    dfuncs::Vector
    #* fluxes graph
    topology::SimpleDiGraph
    #* mtk related
    varinfo::NamedTuple
    paraminfo::NamedTuple
    sys::Union{ODESystem,Nothing}

    function HydroElement(
        name::Symbol;
        funcs::Vector,
        dfuncs::Vector=SimpleFlux[],
        mtk::Bool=true
    )
        topology = build_compute_topology(funcs)
        if !mtk # length(dfuncs) == 0 | 
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

const SurfaceType = :surface
const SoilType = :soil
const FreeWaterType = :freewater

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

function setup_input(
    ele::HydroElement{true};
    input::NamedTuple,
    time::AbstractVector
)
    #* 构建data的插值系统
    itp_eqs = Equation[getproperty(ele.sys, key) ~ @itpfn(key, input[key], time) for key in keys(input)]
    compose_sys = compose(ODESystem(itp_eqs, t; name=Symbol(ele.name, :comp_sys)), ele.sys)
    sys = structural_simplify(compose_sys)
    build_u0 = Pair[]
    for func in filter(func -> func isa AbstractNeuralFlux, ele.funcs)
        func_nn_sys = getproperty(ele.sys, func.param_names)
        push!(build_u0, getproperty(getproperty(func_nn_sys, :input), :u) => zeros(length(get_input_names(func))))
    end
    prob = ODEProblem(sys, build_u0, (time[1], time[end]), [], warn_initialize_determined=true)
    prob
end

function get_sol_0(ele::HydroElement{true};
    input::NamedTuple,
    params::ComponentVector,
    init_states::ComponentVector
)
    ele_input_names = setdiff(get_input_names(ele.funcs),keys(init_states))
    input0 = namedtuple(ele_input_names, [input[nm][1] for nm in ele_input_names])
    init_states_ntp = namedtuple(keys(init_states), [init_states[nm] for nm in keys(init_states)])
    sol_0 = merge(input0, init_states_ntp)
    for func in ele.funcs
        sol_0 = merge(sol_0, func(sol_0, params))
    end
    sol_0
end

function get_mtk_initstates(
    ele::HydroElement{true},
    prob::ODEProblem;
    params::ComponentVector,
    init_states::ComponentVector,
    kw...
)
    #* setup init states
    u0 = [getproperty(ele.sys, nm) => init_states[nm] for nm in keys(init_states) if nm in get_state_names(ele)]
    for func in filter(func -> func isa AbstractNeuralFlux, ele.funcs)
        sol_0 = get_sol_0(ele, input=kw[:input], params=params, init_states=init_states)
        func_nn_sys = getproperty(ele.sys, func.param_names)
        u0 = vcat(u0, [getproperty(getproperty(func_nn_sys, :input), :u)[idx] => sol_0[nm] for (idx, nm) in enumerate(get_input_names(func))])
        # u0 = vcat(u0, [getproperty(getproperty(func_nn_sys, :input), :u) => [sol_0[nm] for nm in get_input_names(func)]])
    end
    u0
end

function get_mtk_params(
    ele::HydroElement{true},
    prob::ODEProblem;
    params::ComponentVector,
    kw...
)
    #* setup init states
    p = Pair[]
    for nm in ModelingToolkit.parameters(ele.sys)
        if contains(string(nm), "₊")
            tmp_nn = split(string(nm), "₊")[1]
            push!(p, getproperty(getproperty(ele.sys, Symbol(tmp_nn)), :p) => Vector(params[Symbol(tmp_nn)]))
        else
            push!(p, getproperty(ele.sys, Symbol(nm)) => params[Symbol(nm)])
        end
    end
    p
end

function setup_prob(
    ele::HydroElement{true},
    prob::ODEProblem;
    params::ComponentVector,
    init_states::ComponentVector,
    kw...
)
    u0 = get_mtk_initstates(ele, prob; params=params, init_states=init_states, kw...)
    p = get_mtk_params(ele, prob; params=params, kw...)
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
    prob = setup_input(ele, input=input[ele_input_names], time=input[:time])
    new_prob = setup_prob(ele, prob, input=input, params=params, init_states=init_states)
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

export HydroElement, add_inputflux!, add_outputflux!, setup_input, solve_prob