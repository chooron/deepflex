"""
HydroElement
"""
struct HydroElement <: AbstractElement
    #* component名称
    name::Symbol
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
        input_names = get_input_names(funcs, dfuncs)
        output_names = get_output_names(funcs)
        state_names = get_state_names(dfuncs)
        param_names = get_param_names(vcat(funcs, dfuncs))

        varinfo, paraminfo = init_var_param(input_names, output_names, state_names, param_names)
        sys = build_ele_system(funcs, dfuncs, varinfo, paraminfo, name=name)
        return new(
            name,
            funcs,
            dfuncs,
            varinfo,
            paraminfo,
            sys
        )
    end
end

function get_ele_io_names(elements::Vector{<:AbstractElement})
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    for ele in elements
        union!(input_names, setdiff(get_input_names(ele.funcs, ele.dfuncs), output_names))
        union!(output_names, get_output_names(ele.funcs))
    end
    input_names, output_names
end

function get_ele_param_names(elements::Vector{<:AbstractElement})
    param_names = Vector{Symbol}()
    for ele in elements
        union!(param_names, get_param_names(vcat(ele.funcs, ele.dfuncs)))
    end
    param_names
end

function get_ele_state_names(elements::Vector{<:AbstractElement})
    state_names = Vector{Symbol}()
    for ele in elements
        union!(state_names, get_state_names(ele.dfuncs))
    end
    state_names
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
        ele_input_names = get_input_names(ele.funcs, ele.dfuncs)
        prob = setup_input(ele, input=fluxes[ele_input_names], time=input[:time])
        new_prob = setup_prob(ele, prob, input=input, params=params, init_states=init_states)
        solved_states = solver(new_prob, get_state_names(ele.dfuncs))
        fluxes = merge(input, solved_states)
    end
    for func in ele.funcs
        fluxes = merge(fluxes, func(fluxes, params))
    end
    fluxes
end

function (elements::AbstractVector{HydroElement})(
    input::NamedTuple,
    pas::ComponentVector;
    solver::AbstractSolver=ODESolver()
)
    fluxes = input
    for ele in elements
        fluxes = merge(fluxes, ele(fluxes, pas, solver=solver))
    end
    fluxes
end

function get_sol_u0(ele::HydroElement, input0::NamedTuple, params::ComponentVector, init_states::NamedTuple)
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
        push!(build_u0, getproperty(getproperty(func_nn_sys, :input), :u) => zeros(eltype(eltype(input)), get_input_names(func)))
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
    ele_input_names = get_input_names(ele.funcs, ele.dfuncs)
    #* 将设定的参数赋予给问题
    #* setup init states
    u0 = Pair[getproperty(ele.sys, nm) => init_states[nm] for nm in keys(init_states) if nm in get_state_names(ele.dfuncs)]
    sol_0 = get_sol_u0(ele, namedtuple(ele_input_names, [input[nm][1] for nm in ele_input_names]),
        params, namedtuple(keys(init_states), [init_states[nm] for nm in keys(init_states)]))
    for func in filter(func -> func isa AbstractNNFlux, ele.funcs)
        func_nn_sys = getproperty(ele.sys, func.param_names)
        u0 = vcat(u0, [getproperty(getproperty(func_nn_sys, :input), :u)[idx] => sol_0[nm] for (idx, nm) in enumerate(get_input_names(func))])
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