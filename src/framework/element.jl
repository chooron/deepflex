
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
    funcs::Vector{AbstractFlux}
    dfuncs::Vector{AbstractFlux}
    lfuncs::Vector{AbstractFlux}
end

function HydroElement(
    ; name::Symbol,
    funcs::Vector{F},
    dfuncs::Vector{F}=SimpleFlux[],
    lfuncs::Vector{LagFlux}=LagFlux[],
) where {F<:AbstractFlux}
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
    name_info = NameInfo(name, input_names, output_names, state_names, param_names)

    return HydroElement(
        name_info,
        funcs,
        dfuncs,
        lfuncs
    )
end

# Element Methods
function get_element_infos(elements::Vector{E}) where {E<:AbstractElement}
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    state_names = Vector{Symbol}()
    param_names = Vector{Symbol}()

    for ele in elements
        union!(input_names, setdiff(ele.name_info.input_names, output_names))
        union!(output_names, ele.output_names)
        union!(param_names, ele.param_names)
        union!(state_names, ele.state_names)
    end
    input_names, output_names, param_names, state_names
end

function solve_prob(
    ele::HydroElement;
    input::NamedTuple,
    params::Union{NamedTuple,ComponentVector},
    init_states::Union{NamedTuple,ComponentVector}
)
    # fit interpolation functions
    ele_input_names = ele.nameinfo.input_names
    ele_state_names = ele.nameinfo.state_names
    itp_dict = namedtuple(ele_input_names, [LinearInterpolation(input[nm], input[:time], extrapolate=true) for nm in ele_input_names])
    function singel_ele_ode_func!(du, u, p, t)
        tmp_input = namedtuple(ele_input_names, [itp_dict[nm](t) for nm in ele_input_names])
        tmp_states = namedtuple(ele_state_names, u)
        tmp_fluxes = merge(tmp_input, tmp_states)
        for func in ele.funcs
            tmp_fluxes = merge(tmp_fluxes, func(tmp_fluxes, params))
        end 
        du[:] = [dfunc(tmp_fluxes, params)[dfunc.output_names] for dfunc in ele.dfuncs]
    end
    prob = ODEProblem(singel_ele_ode_func!, collect(init_states[ele_state_names]), (input[:time][1], input[:time][end]))
    sol = solve(prob, Tsit5(), saveat=input[:time])
    solved_u = hcat(sol.u...)
    sol
    # state_names = collect(keys(init_states))
    # namedtuple(state_names, [solved_u[idx, :] for idx in 1:length(state_names)])
end

function (ele::HydroElement)(
    input::NamedTuple,
    params::Union{NamedTuple,ComponentVector},
    init_states::Union{NamedTuple,ComponentVector}
)
    fluxes = copy(input)
    #* 当存在lagflux时需要提前进行计算
    if length(ele.lfluxes) > 0
        for lf in ele.lfluxes
            fluxes = merge(fluxes, lf(input[lf.input_names]))
        end
    end
    #* 当存在dflux时需要构建ode问题进行计算
    if length(ele.dfuncs) > 0
        solved_states = solve_prob(ele, input=fluxes, params=params,
            init_states=init_states[ele.state_names])
        fluxes = merge(fluxes, solved_states)
    end
    #* 最后就是直接通过flux进行计算
    for func in ele.funcs
        fluxes = merge(fluxes, func(fluxes, params))
    end
    return fluxes
end

function data_itp(t, time::AbstractVector, value::AbstractVector)
    itp = LinearInterpolation(value, time, extrapolate=true)
    itp(t)
end

function build_ele_system(funcs, dfuncs, var_info, param_info; name::Symbol)
    eqs = []
    for func in funcs
        tmp_input = namedtuple(get_input_names(func), [var_info[nm] for nm in get_input_names(func)])
        tmp_param = namedtuple(func.param_names, [param_info[nm] for nm in func.param_names])
        push!(eqs, var_info[first(get_output_names(func))] ~ func.func(tmp_input, tmp_param, func.step_func))
    end
    for dfunc in dfuncs
        tmp_input = namedtuple(get_input_names(dfunc), [var_info[nm] for nm in get_input_names(dfunc)])
        tmp_param = namedtuple(dfunc.param_names, [param_info[nm] for nm in dfunc.param_names])
        push!(eqs, D(var_info[first(get_output_names(dfunc))]) ~ dfunc.func(tmp_input, tmp_param, dfunc.step_func))
    end
    ODESystem(eqs, t; name=name)
end

function build_itp_system(input::NamedTuple, var_info::NamedTuple)
    eqs = []
    for nm in keys(input)
        if nm != :time
            tmp_itp = data_itp(t, input[:time], input[nm])
            func_nm = Symbol(nm, "_itp")
            eval(:($(func_nm)(t) = $(tmp_itp)))
            push!(eqs, eval(Expr(:call, :~, :($(var_info[nm])), :($(func_nm)(t)))))
        end
        #* eval(:(@register_symbolic $(func_nm)(t)))
    end
    ODESystem(eqs, t; name=:itp_sys)
end

function combine_system(ele::HydroElement, system::ODESystem, ipt_system::ODESystem)
    ele_input_names = ele.nameinfo.input_names
    eqs = [eval(Expr(:call, :~, getproperty(system, nm), getproperty(ipt_system, nm))) for nm in ele_input_names]
    compose_sys = compose(ODESystem(eqs, t; name=Symbol(ele.nameinfo.name, :_conn_sys)), system, ipt_system)
    structural_simplify(compose_sys)
end

function init_var_param(nameinfo::NameInfo)
    var_names = vcat(nameinfo.input_names, nameinfo.output_names, nameinfo.state_names)
    var_info = namedtuple(var_names, [first(eval(:(@variables $nm(t)))) for nm in var_names])
    param_info = namedtuple(nameinfo.param_names, [first(eval(:(@parameters $nm))) for nm in nameinfo.param_names])
    var_info, param_info
end

"""
这里只是对funcs和dfuncs两个类型进行求解
"""
function solve_probv2(
    ele::HydroElement;
    input::NamedTuple,
    params::NamedTuple,
    init_states::NamedTuple,
)
    var_info, param_info = init_var_param(ele.nameinfo)
    #* setup data
    itp_sys = build_itp_system(input, var_info)
    base_sys = build_ele_system(ele.funcs, ele.dfuncs, var_info, param_info, name=Symbol(ele.nameinfo.name, :_base_sys))
    sys = combine_system(ele, base_sys, itp_sys)
    #* setup init states
    x0 = Pair{Num,eltype(init_states)}[eval(Expr(:call, :(=>), getproperty(base_sys, nm), init_states[nm])) for nm in ele.nameinfo.state_names]
    #* setup parameters
    p = Pair{Num,eltype(params)}[eval(Expr(:call, :(=>), getproperty(base_sys, Symbol(nm)), params[Symbol(nm)])) for nm in ModelingToolkit.parameters(base_sys)]
    prob = ODEProblem{true,SciMLBase.FullSpecialize}(sys, x0, (1.0, length(input[:time])), p)
    sol = solve(prob, Tsit5(), saveat=input[:time])
    sol
end