struct HydroNode <: AbstractComponent
    name::Symbol
    #* 输入输出名称
    input_names::NamedTuple
    output_names::NamedTuple
    #* 中间状态名称
    state_names::NamedTuple
    #* 参数名称
    param_names::NamedTuple
    #* 单元
    units::AbstractVector{HydroUnit}
    #* 演进模块
    routes::NamedTuple
    #* 节点对应的面积
    area::Number

    function HydroNode(name; unit::AbstractVector{<:AbstractElement}, route::AbstractElement, area::Number=100)
        HydroNode(
            name,
            units=namedtuple([name], [unit]),
            routes=namedtuple([name], [route]),
            area=area
        )
    end

    function HydroNode(name; units::NamedTuple, routes::NamedTuple, area::Number=100)
        param_tuple_list = []
        state_tuple_list = []
        input_names_list = []
        output_names_list = []

        for nm in keys(units)
            ele_list, route = units[nm], routes[nm]
            ele_name_infos = get_element_infos(ele_list)
            route_name_infos = get_element_infos([route])
            push!(param_tuple_list, (unit=ele_name_infos[3], route=route_name_infos[3]))
            push!(state_tuple_list, (unit=ele_name_infos[4], route=route_name_infos[4]))
            push!(input_names_list, ele_name_infos[1])
            push!(output_names_list, unique(vcat(ele_name_infos[2], route_name_infos[2])))
        end

        new(
            name,
            namedtuple(keys(units), input_names_list),
            namedtuple(keys(units), output_names_list),
            namedtuple(keys(units), state_tuple_list),
            namedtuple(keys(units), param_tuple_list),
            units,
            routes,
            area
        )
    end
end

function get_unit_infos(units::AbstractVector{HydroUnit})
    input_names, output_names, state_names, param_names = [], [], [], []
    for unit in units
        push!(input_names, namedtuple([Symbol(unit.name,)], [unit.input_names]))
        push!(output_names, namedtuple([Symbol(unit.name,)], [unit.output_names]))
        push!(state_names, namedtuple([Symbol(unit.name, :_unit)], [unit.state_names]))
        push!(param_names, namedtuple([Symbol(unit.name, :_unit)], [unit.param_names]))
    end
    merge(input_names...), merge(output_names...), merge(param_names...), merge(state_names...)
end

# todo parallel computing
function (node::HydroNode)(
    input::NamedTuple,
    pas::ComponentVector;
    step::Bool=false,
    solver::AbstractSolver=ODESolver()
)
    output_list = []
    for unit in node.units
        tmp_unit_pas = pas[Symbol(unit.name)]
        tmp_ouput = unit(
            input[Symbol(unit.name)],
            tmp_unit_pas[:unit],
            step=step, solver=solver)
        route_output = node.routes[unit.name](tmp_ouput, tmp_unit_pas[:route])
        route_weight_output = route_output[first(node.routes[unit.name].output_names)] .* tmp_unit_pas[:weight]
        push!(output_list, route_weight_output)
    end
    (flow=sum(output_list),)
end