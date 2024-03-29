using DataInterpolations: AbstractInterpolation
import DataInterpolations: derivative
using Symbolics
using Symbolics: Num, unwrap, SymbolicUtils

#* copy from https://github.com/SciML/DataInterpolations.jl/blob/master/ext/DataInterpolationsSymbolicsExt.jl
(interp::AbstractInterpolation)(t::Num) = SymbolicUtils.term(interp, unwrap(t))
SymbolicUtils.promote_symtype(t::AbstractInterpolation, _...) = Real
Base.nameof(interp::AbstractInterpolation) = :Interpolation

@variables t
const D = Differential(t)

struct MTKElement <: AbstractElement
    name::Symbol

    #* attribute
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
    param_names::Vector{Symbol}
    state_names::Vector{Symbol}

    #* functions
    funcs::Vector{AbstractFlux}
    dfuncs::Vector{AbstractFlux}
end

function MTKElement(
    ; name::Symbol,
    funcs::Vector{F},
    dfuncs::Vector{F}
) where {F<:AbstractFlux}
    #* combine the info of func and d_func
    input_names1, output_names, param_names1 = get_func_infos(funcs)
    input_names2, state_names, param_names2 = get_dfunc_infos(dfuncs)
    #* 避免一些中间变量混淆为输入要素
    setdiff!(input_names2, output_names)
    #* 合并两种func的输入要素
    input_names = union(input_names1, input_names2)
    #* 删除输入要素的状态要素
    setdiff!(input_names, state_names)
    param_names = union(param_names1, param_names2)

    return MTKElement(
        name,
        input_names,
        output_names,
        param_names,
        state_names,
        funcs,
        dfuncs,
    )
end

function data_itp(t, time::AbstractVector, value::AbstractVector)
    itp = LinearInterpolation(value, time, extrapolate=true)
    itp(t)
end

function build_ele_system(funcs, dfuncs, var_info, param_info; name::Symbol)
    eqs = []
    #* 根据flux生成方程
    for func in funcs
        #* 筛选出SimpleFlux
        if func isa SimpleFlux
            tmp_input = namedtuple(get_input_names(func), [var_info[nm] for nm in get_input_names(func)])
            tmp_param = namedtuple(func.param_names, [param_info[nm] for nm in func.param_names])
            push!(eqs, var_info[first(get_output_names(func))] ~ func.func(tmp_input, tmp_param, func.step_func))
        end
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

function build_lag_system(input::NamedTuple, var_info::NamedTuple)
    # todo 构建lagsytem前需要根据数据和参数对其进行求解，但是一开始数据未知呀
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

function combine_system(ele::MTKElement, system::ODESystem, ipt_system::ODESystem)
    eqs = [eval(Expr(:call, :~, getproperty(system, nm), getproperty(ipt_system, nm))) for nm in ele.input_names]
    compose_sys = compose(ODESystem(eqs, t; name=Symbol(ele.name, :_conn_sys)), system, ipt_system)
    structural_simplify(compose_sys)
end

function init_var_param(input_names, output_names, param_names, state_names)
    var_names = vcat(input_names, output_names, state_names)
    var_info = namedtuple(var_names, [first(eval(:(@variables $nm(t)))) for nm in var_names])
    param_info = namedtuple(param_names, [first(eval(:(@parameters $nm))) for nm in param_names])
    var_info, param_info
end

function solve_prob(
    ele::MTKElement;
    input::NamedTuple,
    params::NamedTuple,
    init_states::NamedTuple,
)
    var_info, param_info = init_var_param(ele.input_names, ele.output_names, ele.param_names, ele.state_names)
    #* setup data
    base_sys = build_ele_system(ele.funcs, ele.dfuncs, var_info, param_info, name=Symbol(ele.name, :_base_sys))
    itp_sys = build_itp_system(input, var_info)
    #* setup data for lag function
    # lag_sys = build_lag_system(ele.funcs, var_info, input, params)
    sys = combine_system(ele, base_sys, itp_sys)
    #* setup init states
    x0 = Pair{Num,eltype(init_states)}[eval(Expr(:call, :(=>), getproperty(base_sys, nm), init_states[nm])) for nm in ele.state_names]
    #* setup parameters
    p = Pair{Num,eltype(params)}[eval(Expr(:call, :(=>), getproperty(base_sys, Symbol(nm)), params[Symbol(nm)])) for nm in ModelingToolkit.parameters(base_sys)]
    prob = ODEProblem{true,SciMLBase.FullSpecialize}(sys, x0, (1.0, length(input[:time])), p)
    sol = solve(prob, Tsit5(), saveat=input[:time])
    sol
end

function (ele::MTKElement)(
    input::NamedTuple,
    states::NamedTuple,
    params::NamedTuple
)
    fluxes = merge(input, states)
    for func in ele.funcs
        fluxes = merge(fluxes, func(fluxes, params))
    end
    return fluxes
end
