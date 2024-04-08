# todo 像Sciml那样用Unit{false}和Unit{true}来区别基于mtk和非mtk的
struct HydroUnit{E} <: AbstractUnit where {E<:AbstractElement}
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
    elements::AbstractVector{E}
    #* 定义sys
    sys::ODESystem
end

function HydroUnit(name::Symbol;
    elements::Vector{E}
) where {E<:AbstractElement}
    #* 先从element中获取基础信息主要是输出,输入,参数名称等
    input_names, output_names, param_names, state_names = get_element_infos(elements)
    #* 将信息整合到一个类中，便于管理
    nameinfo = NameInfo(name, input_names, output_names, state_names, param_names)
    #* 将这些信息定义成mtk.jl的变量
    sys = build_unit_system(elements, nameinfo)
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

function get_element_infos(elements::Vector{E}) where {E<:AbstractElement}
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
    elements::AbstractVector{E},
    nameinfo::NamedTuple,
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
                push!(eqs, getproperty(tmp_ele1.base_sys, nm) ~ getproperty(tmp_ele2.base_sys, nm))
            end
        end
    end
    compose(ODESystem(eqs, t; name=Symbol(nameinfo.name, :_sys)), [ele.sys for ele in elements]...)
end

function setup_input(unit::HydroUnit; input::NamedTuple, time::AbstractVector)
    #* 首先构建data的插值系统
    eqs = []
    itp_sys = build_itp_system(input[ele.input_names], time, ele.varinfo, name=ele.name)
    for ele in unit.elements
        for nm in filter(nm -> nm in ele.input_names, keys(input))
            push!(eqs, getproperty(ele.sys, nm) ~ getproperty(itp_sys, nm))
        end
    end
    compose_sys = compose(ODESystem(eqs, t; name=Symbol(ele.name, :comp_sys)), unit.sys, itp_sys)
    structural_simplify(compose_sys)
end


function (unit::HydroUnit)(
    input::NamedTuple,
    params::NamedTuple,
    init_states::NamedTuple,
    step::Bool=false
)
    if step
        return _step_forward(unit, input, params, init_states)
    else
        return _whole_forward(unit, input, params, init_states)
    end
end

function _step_forward(
    unit::HydroUnit,
    input::NamedTuple,
    params::NamedTuple,
    init_states::NamedTuple,
)
    # * This function is calculated element by element
    fluxes = input
    for ele in unit.elements
        fluxes = merge(fluxes, ele(fluxes, params, init_states))
    end
    fluxes
end

function _whole_forward(
    unit::HydroUnit,
    input::NamedTuple,
    params::NamedTuple,
    init_states::NamedTuple,
)
    sys = setup_input(unit, input=input[unit.input_names], time=input[:time])
    solved_states = solve_prob(unit.elements,
        sys=sys, input=fluxes, params=params,
        init_states=init_states[ele.state_names]
    )
    solved_states
end

function get_all_luxnnflux(unit::HydroUnit)
    luxnn_tuple = namedtuple()
    for ele in unit.elements
        merge!(luxnn_tuple, get_all_luxnnflux(ele))
    end
    luxnn_tuple
end