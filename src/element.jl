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
    lfuncs::Vector
    #* mtk related
    varinfo::NamedTuple
    paraminfo::NamedTuple
    sys::Union{ODESystem,Nothing}
end

function HydroElement(
    ; name::Symbol,
    funcs::Vector,
    dfuncs::Vector=[],
    lfuncs::Vector=[],
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
    varinfo, paraminfo = init_var_param(nameinfo)
    sys = build_ele_system(funcs, dfuncs, varinfo, paraminfo, name=name)
    return HydroElement(
        name,
        input_names,
        output_names,
        state_names,
        param_names,
        funcs,
        dfuncs,
        lfuncs,
        varinfo,
        paraminfo,
        sys
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

function (ele::HydroElement)(
    input::NamedTuple,
    params::Union{NamedTuple,ComponentVector},
    init_states::Union{NamedTuple,ComponentVector},
    solver::AbstractSolver
)
    # todo 这里需要添加input，params，init_states的参数校核
    fluxes = input
    #* 当存在lagflux时需要提前进行计算
    if length(ele.lfuncs) > 0
        for lf in ele.lfuncs
            fluxes = merge(fluxes, lf(input, params))
        end
    end
    #* 当存在dflux时需要构建ode问题进行计算
    if length(ele.dfuncs) > 0
        #* 根据方程构建基础系统
        sys = setup_input(ele, input=input, time=input[:time])
        #* 设定系统参数构建ode problem
        prob = build_prob(ele,
            sys=sys, input=input, params=params,
            init_states=init_states[ele.state_names])
        solved_states = solver(prob, ele.state_names)
        fluxes = merge(fluxes, solved_states)
    end
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
    ele::HydroElement;
    sys::ODESystem,
    input::NamedTuple,
    params::NamedTuple,
    init_states::NamedTuple,
)
    #* setup init states
    x0 = Pair{Num,eltype(init_states)}[getproperty(ele.sys, nm) => init_states[nm] for nm in keys(init_states)]
    #* setup parameters
    p = Pair{Num,eltype(params)}[getproperty(ele.sys, Symbol(nm)) => params[Symbol(nm)] for nm in ModelingToolkit.parameters(ele.sys)]
    #* build problem
    ODEProblem{true,SciMLBase.FullSpecialize}(sys, x0, (input[:time][1], input[:time][end]), p)
end

function build_prob(
    elements::AbstractVector{HydroElement};
    sys::ODESystem,
    input::NamedTuple,
    params::NamedTuple,
    init_states::NamedTuple,
)
    #* setup init states and parameters
    x0 = Pair{Num,eltype(init_states)}[]
    p = Pair{Num,eltype(params)}[]
    for ele in elements
        for nm in filter(nm -> nm in ele.state_names, keys(init_states))
            push!(x0, getproperty(ele.sys, Symbol(nm)) => init_states[Symbol(nm)])
        end
        for nm in ModelingToolkit.parameters(ele.sys)
            push!(p, getproperty(ele.sys, Symbol(nm)) => params[Symbol(nm)])
        end
    end
    #* prepare problem args
    ODEProblem{true,SciMLBase.FullSpecialize}(sys, x0, (input[:time][1], input[:time][end]), p)
end

function build_prob(
    ele::HydroElement;
    input::NamedTuple,
    params::Union{NamedTuple,ComponentVector},
    init_states::Union{NamedTuple,ComponentVector}
)
    # fit interpolation functions
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