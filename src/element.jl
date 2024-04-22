"""
HydroElement
"""
struct HydroElement <: AbstractElement
    #* component名称
    name::Symbol
    #* 输入输出名称
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
    #* 中间状态名称
    state_names::Vector{Symbol}
    #* 参数名称
    param_names::Vector{Symbol}
    #* fluxes and dfluxes
    funcs::Vector
    dfuncs::Vector
    #* mtk related
    varinfo::NamedTuple
    paraminfo::NamedTuple
    sys::Union{ODESystem,Nothing}
    prob::Union{ODEProblem,Nothing}
end

struct RouteElement <: AbstractElement
    #* component名称
    name::Symbol
    #* 输入输出名称
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
    #* 参数名称
    param_names::Vector{Symbol}
    #* fluxes (混有simpleflux和RouteFlux)
    funcs::Vector
end

function HydroElement(
    ; name::Symbol,
    funcs::Vector,
    dfuncs::Vector=SimpleFlux[],
)
    # combine the info of func and d_func
    input_names1, output_names, param_names1 = get_func_infos(funcs)
    input_names2, state_names, param_names2 = get_dfunc_infos(dfuncs)
    # 避免一些中间变量混淆为输入要素
    setdiff!(input_names2, output_names)
    # 合并两种func的输入要素
    input_names = union(input_names1, input_names2)
    # 删除输入要素的状态要素
    setdiff!(input_names, state_names)
    # 合并两种类型函数的参数
    param_names = union(param_names1, param_names2)
    # 条件判断，有时候不需要构建
    varinfo, paraminfo = init_var_param(input_names, output_names, state_names, param_names)
    sys = build_ele_system(funcs, dfuncs, varinfo, paraminfo, name=name)
    prob = nothing
    return HydroElement(
        name,
        input_names,
        output_names,
        state_names,
        param_names,
        funcs,
        dfuncs,
        varinfo,
        paraminfo,
        sys,
        prob
    )
end

function RouteElement(;
    name::Symbol,
    funcs::Vector)

    # funcs including both SimpleFLux and RouteFlux
    input_names, output_names, param_names = get_func_infos(funcs)

    return RouteElement(
        name,
        input_names,
        output_names,
        param_names,
        funcs,
    )
end

# Element Methods
function get_all_luxnnflux(ele::HydroElement)
    luxnn_tuple = namedtuple()
    for func in vcat(ele.funcs, ele.dfuncs)
        if func isa AbstractNNFlux
            merge!(luxnn_tuple, namedtuple([func.param_names], [func]))
        end
    end
    luxnn_tuple
end


# 求解并计算
function (ele::HydroElement)(
    input::NamedTuple,
    pas::ComponentVector;
    solver::AbstractSolver=ODESolver()
)
    params, init_states = pas[:params], pas[:initstates]
    fluxes = input
    if length(ele.dfuncs) > 0
        solved_states = solve_prob(ele, input, params, init_states, solver=solver)
        fluxes = merge(input, solved_states)
    end
    for func in ele.funcs
        fluxes = merge(fluxes, func(fluxes, params))
    end
    fluxes
end

function (ele::RouteElement)(;
    input::NamedTuple,
    pas::ComponentVector,
)
    params = pas[:params]
    params = params isa AbstractVector ? ComponentVector() : params
    fluxes = input
    for func in ele.funcs
        fluxes = merge(fluxes, func(fluxes, params))
    end
    fluxes
end

function (ele::RouteElement)(
    input::NamedTuple,
    params::ComponentVector,
)
    # todo 这里需要添加input，params，init_states的参数校核
    fluxes = input
    #* 最后就是直接通过flux进行计算
    for func in ele.funcs
        fluxes = merge(fluxes, func(fluxes, params))
    end
    return fluxes
end

function setup_input(
    ele::HydroElement;
    input::NamedTuple,
    time::AbstractVector
)
    #* 首先构建data的插值系统
    itp_sys = build_itp_system(input[ele.input_names], time, ele.varinfo, name=ele.name)
    #* 连接系统
    eqs = [getproperty(ele.sys, nm) ~ getproperty(itp_sys, nm) for nm in ele.input_names]
    compose_sys = compose(ODESystem(eqs, t; name=Symbol(ele.name, :comp_sys)), ele.sys, itp_sys)
    structural_simplify(compose_sys)
end

function build_prob(
    ele::HydroElement,
    sys::ODESystem;
    input::NamedTuple,
    params::ComponentVector,
    init_states::ComponentVector,
)
    #* setup init states
    x0 = Pair{Num,eltype(init_states)}[getproperty(ele.sys, nm) => init_states[nm] for nm in keys(init_states)]
    #* setup parameters
    p = Pair{Num,eltype(params)}[getproperty(ele.sys, Symbol(nm)) => params[Symbol(nm)] for nm in ModelingToolkit.parameters(ele.sys)]
    #* build problem
    ODEProblem{true,SciMLBase.FullSpecialize}(sys, x0, (input[:time][1], input[:time][end]), p)
end


function build_prob(
    ele::HydroElement;
    input::NamedTuple,
    params::Union{NamedTuple,ComponentVector},
    init_states::Union{NamedTuple,ComponentVector}
)
    itp_dict = namedtuple(ele.input_names, [LinearInterpolation(input[nm], input[:time], extrapolate=true) for nm in ele.input_names])
    function singel_ele_ode_func!(du, u, p, t)
        tmp_input = namedtuple(ele.input_names, [itp_dict[nm](t) for nm in ele.input_names])
        tmp_states = namedtuple(ele.state_names, u)
        tmp_fluxes = merge(tmp_input, tmp_states)
        for func in ele.funcs
            tmp_fluxes = merge(tmp_fluxes, func(tmp_fluxes, params))
        end
        du[:] = [dfunc(tmp_fluxes, params)[dfunc.output_names] for dfunc in ele.dfuncs]
    end
    prob = ODEProblem(singel_ele_ode_func!, collect(init_states[ele.state_names]), (input[:time][1], input[:time][end]))
    prob
end

function solve_prob(
    ele::HydroElement,
    input::NamedTuple,
    params::ComponentVector,
    init_states::ComponentVector;
    solver::AbstractSolver=ODESolver()
)
    #* 当存在dflux时需要构建ode问题进行计算
    #* 根据方程构建基础系统
    sys = setup_input(ele, input=input, time=input[:time])
    #* 设定系统参数构建ode problem
    prob = build_prob(ele, sys,
        input=input, params=params,
        init_states=init_states[ele.state_names])
    solver(prob, ele.state_names)
end