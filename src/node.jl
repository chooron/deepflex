struct HydroNode <: AbstractComponent
    name::Symbol
    #* 单元
    units::NamedTuple
    #* 节点对应的面积
    area::Number

    function HydroNode(name; units::Union{NamedTuple, AbstractVector{<:AbstractElement}}, area::Number=100)
        units = units isa AbstractVector ? namedtuple([name], [units]) : units
        new(
            name,
            units,
            area
        )
    end
end

function get_ele_io_names(units::NamedTuple)
    input_names_list = []
    output_names_list = []
    for unit in keys(units)
        ele_input_names, ele_output_names = get_ele_io_names(unit)
        push!(input_names_list, ele_input_names)
        push!(output_names_list, ele_output_names)
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
    input::ComponentVector,
    pas::ComponentVector;
    solver::AbstractSolver=ODESolver()
)
    output_list = []
    for nm in keys(node.units)
        unit = node.units[nm]
        tmp_ouput = unit(input[nm], pas[nm], solver=solver)
        # println(namedtuple(keys(tmp_ouput), [length(tmp_ouput[nm]) for nm in keys(tmp_ouput)]))
        push!(output_list, tmp_ouput[:flow] .* pas[nm][:weight])
    end
    (flow=sum(output_list),)
end

# function setup_input!(node::HydroNode; input::ComponentVector)
#     for ele in node.

# end