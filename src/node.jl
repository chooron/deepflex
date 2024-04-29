struct HydroNode <: AbstractComponent
    name::Symbol
    #* 单元
    units::NamedTuple
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
        new(
            name,
            units,
            routes,
            area
        )
    end
end

function get_ele_io_names(units::NamedTuple, routes::NamedTuple)
    input_names_list = []
    output_names_list = []
    for nm in keys(units)
        ele_input_names, ele_output_names = get_ele_io_names(units[nm])
        push!(input_names_list, ele_input_names)
        push!(output_names_list, unique(vcat(ele_output_names, get_ele_io_names([routes[nm]])[2])))
    end
    namedtuple(keys(units), input_names_list), namedtuple(keys(units), output_names_list)
end

function get_ele_param_names(units::NamedTuple)
    param_tuple_list = []
    for nm in keys(units)
        push!(param_tuple_list, (unit=get_param_names(units[nm]), route=get_param_names([route])))
    end
    param_tuple_list
end

function get_ele_state_names(units::NamedTuple)
    state_tuple_list = []
    for nm in keys(units)
        push!(state_tuple_list, (unit=get_state_names(units[nm]), route=get_state_names([route])))
    end
    state_tuple_list
end


# todo parallel computing
function (node::HydroNode)(
    input::NamedTuple,
    pas::ComponentVector;
    solver::AbstractSolver=ODESolver()
)
    output_list = []
    for nm in keys(node.units)
        unit = node.units[nm]
        tmp_unit_pas = pas[nm]
        tmp_ouput = unit(
            input[nm],
            tmp_unit_pas[:unit],
            solver=solver)
        route_output = node.routes[nm](tmp_ouput, tmp_unit_pas[:route])
        #todo 这里需要改一改
        route_weight_output = route_output[first(get_ele_io_names([node.routes[nm]])[2])] .* tmp_unit_pas[:weight]
        push!(output_list, route_weight_output)
    end
    (flow=sum(output_list),)
end