struct HydroNode <: AbstractNode
    name::Symbol
    #* 子单元名称
    subnames::Vector{Symbol}
    #* 子单元 NamedTuple{Symbol,{HydroUnit}}
    units::NamedTuple
    #* RunoffRouting Module NamedTuple{Symbol,HydroElement}
    routes::NamedTuple
    #* 节点对应的面积
    area::Number
    
    function HydroNode(name; units::Vector{<:HydroUnit}, routes::Vector{<:HydroElement}, area::Number=100.0)
        subnames = [unit.name for unit in units]
        units_ntp = namedtuple([unit.name for unit in units], units)
        routes_ntp = namedtuple([route.name for route in routes], routes)

        @assert keys(units_ntp) == keys(routes_ntp)
        new(
            name,
            subnames,
            units_ntp,
            routes_ntp,
            area
        )
    end
end

function (node::HydroNode)(
    input::NamedTuple,
    pas::ComponentVector;
    solver::AbstractSolver=ODESolver()
)
    output_list = []

    @threads for idx in 1:length(node.subnames)
        subname = node.subnames[idx]
        tmp_input = input[subname]
        tmp_unit_ouput = node.units[idx](tmp_input, pas[subname], solver=solver)
        tmp_route_output = node.routes[idx](tmp_unit_ouput, pas[subname])
        push!(output_list, tmp_route_output[:flow] .* pas[subname][:weight])
    end
    (flow=sum(output_list),)
end

export HydroNode