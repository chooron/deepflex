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
end

function HydroNode(name; unit::HydroUnit, route::AbstractElement)
    HydroNode(
        name,
        units=[unit],
        routes=namedtuple([name], [route]),
    )
end

function HydroNode(name; units::AbstractVector{HydroUnit}, routes::NamedTuple, area::Number=100)
    input_names, output_names, param_names, state_names = get_unit_infos(units)

    #* 合并参数名称
    param_names = merge(param_names, (weights=[unit.name for unit in units],))
    for (nm, route) in pairs(routes)
        param_names = merge(param_names, namedtuple([Symbol(nm, :_route)], [route.param_names]))
    end

    HydroNode(
        name,
        input_names,
        output_names,
        state_names,
        param_names,
        units,
        routes,
        area
    )
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