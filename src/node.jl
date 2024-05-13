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
    
    function HydroNode(name; units::Vector{HydroUnit}, routes::Vector{HydroElement}, area::Number=100.0)
        units_ntp = namedtuple([unit.name for unit in units], units)
        routes_ntp = namedtuple([route.name for route in routes], routes)
        @assert keys(units_ntp) == keys(routes_ntp)
        subnames = collect(keys(units))

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
    @threads for (idx, nm) in enumerate(node.subnames)
        tmp_input = input[nm]
        tmp_unit_ouput = node.units[idx](tmp_input, pas[nm], solver=solver)
        tmp_route_output = node.routes[idx](tmp_unit_ouput, pas[nm])
        push!(output_list, tmp_route_output[:flow] .* pas[nm][:weight])
    end
    (flow=sum(output_list),)
end

export HydroNode