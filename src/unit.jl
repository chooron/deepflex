# todo 像Sciml那样用Unit{false}和Unit{true}来区别基于mtk和非mtk的
struct HydroUnit <: AbstractUnit
    #* component名称
    name::Symbol
    #* 输入输出名称
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
    #* 中间状态名称
    state_names::Vector{Symbol}
    #* 参数名称
    param_names::Vector{Symbol}
    #* 产流计算元素组
    elements::AbstractVector
    #* 定义sys
    sys::ODESystem
end

function HydroUnit(
    name::Symbol;
    elements::Vector{E}
) where {E<:AbstractElement}
    #* 先从element中获取基础信息主要是输出,输入,参数名称等
    input_names, output_names, param_names, state_names = get_element_infos(elements)
    #* 将这些信息定义成mtk.jl的变量
    sys = build_unit_system(elements, name=name)
    HydroUnit(
        name,
        input_names,
        output_names,
        state_names,
        param_names,
        elements,
        sys
    )
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

function build_unit_system(
    elements::AbstractVector{E};
    name::Symbol,
) where {E<:AbstractElement}
    eqs = []
    elements = elements
    #* 这里是遍历两次这个elements，就是通过对比两两数组的共同flux名称，然后将名称相同的名称构建方程
    for tmp_ele1 in elements
        #* 连接element之间的变量
        for tmp_ele2 in filter(ele -> ele.name != tmp_ele1.name, elements)
            share_var_names = intersect(
                vcat(tmp_ele1.input_names, tmp_ele1.output_names, tmp_ele1.state_names),
                vcat(tmp_ele2.input_names, tmp_ele2.output_names, tmp_ele2.state_names)
            )
            for nm in share_var_names
                push!(eqs, getproperty(tmp_ele1.sys, nm) ~ getproperty(tmp_ele2.sys, nm))
            end
        end
    end
    compose(ODESystem(eqs, t; name=Symbol(name, :_sys)), [ele.sys for ele in elements]...)
end

function setup_input(
    unit::HydroUnit;
    input::NamedTuple,
    time::AbstractVector
)
    #* 首先构建data的插值系统
    eqs = Equation[]
    unit_varinfo = merge([ele.varinfo for ele in unit.elements]...)
    itp_sys = build_itp_system(input[unit.input_names], time, unit_varinfo, name=unit.name)
    for ele in unit.elements
        for nm in filter(nm -> nm in ele.input_names, keys(input))
            push!(eqs, getproperty(getproperty(unit.sys, ele.sys.name), nm) ~ getproperty(itp_sys, nm))
        end
    end
    compose_sys = compose(ODESystem(eqs, t; name=Symbol(unit.name, :_comp_sys)), unit.sys, itp_sys)
    structural_simplify(compose_sys)
end

function build_prob(
    unit::HydroUnit,
    sys::ODESystem;
    input::NamedTuple,
    params::ComponentVector,
    init_states::ComponentVector,
)
    #* setup init states and parameters
    x0 = Pair{Num,eltype(init_states)}[]
    p = Pair{Num,eltype(params)}[]

    for ele in unit.elements
        for nm in filter(nm -> nm in ele.state_names, keys(init_states))
            push!(x0, getproperty(getproperty(unit.sys, ele.sys.name), Symbol(nm)) => init_states[Symbol(nm)])
        end
        for nm in ModelingToolkit.parameters(ele.sys)
            push!(p, getproperty(getproperty(unit.sys, ele.sys.name), Symbol(nm)) => params[Symbol(nm)])
        end
    end

    #* build problem
    prob = ODEProblem{true,SciMLBase.FullSpecialize}(sys, x0, (input[:time][1], input[:time][end]), p)
    println(parameter_values(prob))
    
    prob
end

function (unit::HydroUnit)(
    input::NamedTuple,
    pas::ComponentVector;
    step::Bool=false,
    solver::AbstractSolver=ODESolver()
)
    if step
        return _step_forward(unit, input, pas, solver)
    else
        return _whole_forward(unit, input, pas, solver)
    end
end

function _step_forward(
    unit::HydroUnit,
    input::NamedTuple,
    pas::ComponentVector,
    solver::AbstractSolver,
)
    # * This function is calculated element by element
    fluxes = input
    for ele in unit.elements
        fluxes = merge(fluxes, ele(fluxes, pas, solver=solver))
    end
    fluxes
end

function _whole_forward(
    unit::HydroUnit,
    input::NamedTuple,
    pas::ComponentVector,
    solver::AbstractSolver,
)
    params, init_states = pas[:params], pas[:initstates]
    sys = setup_input(unit, input=input[unit.input_names], time=input[:time])
    prob = build_prob(unit, sys,
        input=input, params=params,
        init_states=init_states[unit.state_names]
    )
    solved_states = solver(prob, unit.state_names)
    fluxes = merge(input, solved_states)
    for ele in unit.elements
        for func in ele.funcs
            fluxes = merge(fluxes, func(fluxes, params))
        end
    end
    fluxes
end

function get_all_luxnnflux(unit::HydroUnit)
    luxnn_tuple = namedtuple()
    for ele in unit.elements
        merge!(luxnn_tuple, get_all_luxnnflux(ele))
    end
    luxnn_tuple
end