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

    function HydroElement(
        name::Symbol;
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
        return new(
            name,
            input_names,
            output_names,
            state_names,
            param_names,
            funcs,
            dfuncs,
            varinfo,
            paraminfo,
            sys
        )
    end

end

function get_element_infos(elements::Vector{HydroElement})
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    state_names = Vector{Symbol}()
    param_names = Vector{Symbol}()

    for ele in elements
        union!(input_names, setdiff(ele.input_names, output_names))
        union!(output_names, ele.output_names)
        union!(param_names, ele.param_names)
        union!(state_names, ele.state_names)
    end
    input_names, output_names, param_names, state_names
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
    params = pas[:params] isa Vector ? ComponentVector() : pas[:params]
    fluxes = input
    if length(ele.dfuncs) > 0
        init_states = pas[:initstates] isa Vector ? ComponentVector() : pas[:initstates]
        prob = setup_input(ele, input=fluxes[ele.input_names], time=input[:time])
        new_prob = setup_prob(ele, prob, input=input, params=params, init_states=init_states)
        solved_states = solver(new_prob, ele.state_names)
        fluxes = merge(input, solved_states)
    end
    for func in ele.funcs
        fluxes = merge(fluxes, func(fluxes, params))
    end
    fluxes
end

function get_sol_u0(ele, input0, params, init_states)
    #* 跳过常微分方程求解直接计算各中间变量的初始状态
    input = merge(input0, init_states)
    for func in ele.funcs
        input = merge(input, func(input, params))
    end
    input
end

function setup_input(
    ele::HydroElement;
    input::NamedTuple,
    time::AbstractVector
)
    #* 构建data的插值系统
    itp_eqs = Equation[getproperty(ele.sys, nm) ~ @itpfn(nm, ip, time) for (nm, ip) in pairs(input)]
    compose_sys = compose(ODESystem(itp_eqs, t; name=Symbol(ele.name, :comp_sys)), ele.sys)
    sys = structural_simplify(compose_sys)
    build_u0 = Pair[]
    for func in filter(func -> func isa AbstractNNFlux, ele.funcs)
        func_nn_sys = getproperty(ele.sys, func.param_names)
        push!(build_u0, getproperty(getproperty(func_nn_sys, :input), :u) => zeros(eltype(eltype(input)), length(func.input_names)))
    end
    prob = ODEProblem(sys, build_u0, (time[1], time[end]), [], warn_initialize_determined=true)
    prob
end

function setup_prob(
    ele::HydroElement,
    prob::ODEProblem;
    input::NamedTuple,
    params::ComponentVector,
    init_states::ComponentVector,
)
    #* 将设定的参数赋予给问题
    #* setup init states
    u0 = Pair[getproperty(ele.sys, nm) => init_states[nm] for nm in keys(init_states) if nm in ele.state_names]
    sol_0 = get_sol_u0(ele, namedtuple(ele.input_names, [input[nm][1] for nm in ele.input_names]),
        params, namedtuple(keys(init_states), [init_states[nm] for nm in keys(init_states)]))
    for func in filter(func -> func isa AbstractNNFlux, ele.funcs)
        func_nn_sys = getproperty(ele.sys, func.param_names)
        u0 = vcat(u0, [getproperty(getproperty(func_nn_sys, :input), :u)[idx] => sol_0[nm] for (idx, nm) in enumerate(func.input_names)])
    end
    #* setup parameters
    p = Pair[]
    for nm in ModelingToolkit.parameters(ele.sys)
        if contains(string(nm), "₊")
            tmp_nn = split(string(nm), "₊")[1]
            push!(p, getproperty(getproperty(ele.sys, Symbol(tmp_nn)), :p) => Vector(params[Symbol(tmp_nn)]))
        else
            push!(p, getproperty(ele.sys, Symbol(nm)) => params[Symbol(nm)])
        end
    end
    new_prob = remake(prob, p=p, u0=u0)
    new_prob
end