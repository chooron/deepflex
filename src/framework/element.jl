
struct NameInfo
    #* component名称
    name::Symbol
    #* 输入输出名称
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
    #* 中间状态名称
    state_names::Vector{Symbol}
    #* 参数名称
    param_names::Vector{Symbol}
end

"""
HydroElement
"""
struct HydroElement <: AbstractElement
    #* 名称信息
    nameinfo::NameInfo
    #* fluxes and dfluxes
    funcs::Vector
    dfuncs::Vector
    lfuncs::Vector
    #* mtk related
    varinfo::NamedTuple
    paraminfo::NamedTuple
    base_sys::ODESystem
end

function HydroElement(
    ; name::Symbol,
    funcs::Vector,
    dfuncs::Vector=SimpleFlux[],
    lfuncs::Vector=LagFlux[],
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
    # 检查lag fluxes的输入只能是input name，因为lag flux通常是需要在最开始进行计算
    nameinfo = NameInfo(name, input_names, output_names, state_names, param_names)

    varinfo, paraminfo = init_var_param(nameinfo)
    base_sys = build_ele_system(funcs, dfuncs, nameinfo, varinfo, paraminfo)
    return HydroElement(
        nameinfo,
        funcs,
        dfuncs,
        lfuncs,
        varinfo,
        paraminfo,
        base_sys
    )
end

# Element Methods
function get_element_infos(elements::Vector{E}) where {E<:AbstractElement}
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    state_names = Vector{Symbol}()
    param_names = Vector{Symbol}()

    for ele in elements
        union!(input_names, setdiff(ele.nameinfo.input_names, output_names))
        union!(output_names, ele.nameinfo.output_names)
        union!(param_names, ele.nameinfo.param_names)
        union!(state_names, ele.nameinfo.state_names)
    end
    input_names, output_names, param_names, state_names
end

function get_all_luxnnflux(ele::HydroElement)
    luxnn_tuple = namedtuple()
    for func in vcat(ele.funcs, ele.dfuncs, ele.lfuncs)
        if func isa AbstractNNFlux
            merge!(luxnn_tuple, namedtuple([func.param_names], [func]))
        end
    end
    luxnn_tuple
end

function (ele::HydroElement)(
    input::NamedTuple,
    params::Union{NamedTuple,ComponentVector},
    init_states::Union{NamedTuple,ComponentVector}
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
        # 根据方程构建基础系统
        sys = setup_data(ele, input=input, time=input[:time])
        solved_states = solve_prob(ele,
            sys=sys, input=fluxes, params=params,
            init_states=init_states[ele.nameinfo.state_names]
        )
        fluxes = merge(fluxes, solved_states)
    end
    #* 最后就是直接通过flux进行计算
    for func in ele.funcs
        fluxes = merge(fluxes, func(fluxes, params))
    end
    return fluxes
end

function init_var_param(nameinfo::NameInfo)
    var_names = vcat(nameinfo.input_names, nameinfo.output_names, nameinfo.state_names)
    varinfo = namedtuple(var_names, [first(@variables $nm(t)) for nm in var_names])
    paraminfo = namedtuple(nameinfo.param_names, [first(@parameters $nm) for nm in nameinfo.param_names])
    varinfo, paraminfo
end

function build_ele_system(funcs, dfuncs, nameinfo, varinfo::NamedTuple, paraminfo::NamedTuple)
    eqs = Equation[]
    for func in funcs
        tmp_input = namedtuple(get_input_names(func), [varinfo[nm] for nm in get_input_names(func)])
        tmp_param = namedtuple(get_param_names(func), [paraminfo[nm] for nm in get_param_names(func)])
        push!(eqs, varinfo[first(get_output_names(func))] ~ func(tmp_input, tmp_param)[first(get_output_names(func))])
    end
    for dfunc in dfuncs
        tmp_input = namedtuple(get_input_names(dfunc), [varinfo[nm] for nm in get_input_names(dfunc)])
        tmp_param = namedtuple(get_param_names(dfunc), [paraminfo[nm] for nm in get_param_names(dfunc)])
        push!(eqs, D(varinfo[first(get_output_names(dfunc))]) ~ dfunc(tmp_input, tmp_param)[first(get_output_names(dfunc))])
    end
    ODESystem(eqs, t; name=Symbol(nameinfo.name, :base_sys))
end

function build_itp_system(input::NamedTuple, time::AbstractVector, varinfo::NamedTuple; name::Symbol)
    eqs = Equation[]
    for nm in keys(input)
        tmp_itp = LinearInterpolation(input[nm], time, extrapolate=true)
        push!(eqs, varinfo[nm] ~ tmp_itp(t))
    end
    ODESystem(eqs, t; name=Symbol(name, :itp_sys))
end

function setup_data(ele::HydroElement; input::NamedTuple, time::AbstractVector)
    #* 首先构建data的插值系统
    itp_sys = build_itp_system(input[ele.nameinfo.input_names], time, ele.varinfo, name=ele.nameinfo.name)
    eqs = [getproperty(ele.base_sys, nm) ~ getproperty(itp_sys, nm) for nm in ele.nameinfo.input_names]
    compose_sys = compose(ODESystem(eqs, t; name=Symbol(ele.nameinfo.name, :comp_sys)), ele.base_sys, itp_sys)
    structural_simplify(compose_sys)
end

function solve_prob(
    ele::HydroElement;
    sys::ODESystem,
    input::NamedTuple,
    params::NamedTuple,
    init_states::NamedTuple,
)
    #* setup init states
    x0 = Pair{Num,eltype(init_states)}[getproperty(ele.base_sys, nm) => init_states[nm] for nm in keys(init_states)]
    #* setup parameters
    p = Pair{Num,eltype(params)}[getproperty(ele.base_sys, Symbol(nm)) => params[Symbol(nm)] for nm in ModelingToolkit.parameters(ele.base_sys)]
    prob = ODEProblem{true,SciMLBase.FullSpecialize}(sys, x0, (1.0, length(input[:time])), p)
    sol = solve(prob, Tsit5(), saveat=input[:time])
    sol_u = hcat(sol.u...)
    namedtuple(keys(init_states), [sol_u[i, :] for i in 1:size(sol_u)[1]])
end