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
    routes::NamedTuple
end

function HydroNode(name; unit::HydroUnit, route::AbstractElement)
    HydroNode(
        name,
        units=[unit],
        routes=namedtuple([name], [route]),
    )
end

function HydroNode(name; units::AbstractVector{HydroUnit}, routes::NamedTuple)
    input_names, output_names, param_names, state_names = get_unit_infos(units)

    #* 合并参数名称
    merge!(param_names, (weights=[Symbol(unit.name, :_weight) for unit in units],))
    route_params = namedtuple()
    for (unit, route) in pairs(routes)
        merge!(route_params, namedtuple([Symbol(unit.name, :_route)], route.param_names))
    end
    merge!(param_names, route_params)

    HydroNode(
        name,
        input_names,
        output_names,
        state_names,
        param_names,
        units,
        routes,
    )
end

function get_unit_infos(units::AbstractVector{HydroUnit})
    input_names = namedtuple()
    output_names = namedtuple()
    state_names = namedtuple()
    param_names = namedtuple()
    for unit in units
        merge!(input_names, namedtuple([Symbol(unit.name)], route.input_names))
        merge!(output_names, namedtuple([Symbol(unit.name)], route.output_names))
        merge!(state_names, namedtuple([Symbol(unit.name)], route.state_names))
        merge!(param_names, namedtuple([Symbol(unit.name, :_unit)], route.param_names))
    end
    input_names, output_names, param_names, state_names
end

function inner_routing(flux::AbstractVector, routefunc::AbstractElement, params::NamedTuple, init_states::NamedTuple)
    routefunc(flux, params, init_states)
end

# todo parallel computing
function (node::Node)(input::NamedTuple, params::NamedTuple, init_states::NamedTuple)
    unit_names = [u.name for u in node.units]
    unit_ouputs = namedtuple(unit_names, [
        begin
            tmp_ouput = u(input, params, init_states)
            route_output = node.routes[u.name](tmp_ouput, params, init_states) .* node.weights[u.name]
            route_output
        end
        for u in node.units
    ])
    unit_ouputs
end