"""
$(TYPEDEF)
A basic hydrological calculation node can represent a complete conceptual hydrological model including vertical and horizontal calculations,
 can contain different vertical calculation units, and can calculate distributed and semi-distributed hydrological models in scope.
# Fields
$(FIELDS)
# Example
```
using LumpedHydro
using LumpedHydro.ExpHydro

model = LumpedHydro.HydroNode(
    :exphydro_node,
    units=[LumpedHydro.ExpHydro.Unit(name=:exphydro1, mtk=true, step=false), LumpedHydro.ExpHydro.Unit(name=:exphydro2, mtk=true, step=false)],
    routes=[LumpedHydro.ExpHydro.Route(name=:exphydro1), LumpedHydro.ExpHydro.Route(name=:exphydro2)]
)
```
"""
struct HydroNode <: AbstractNode
    "hydrological computation node name"
    name::Symbol
    "sub-hydrological computation unit names"
    subnames::Vector{Symbol}
    "sub-hydrological computation unit: (unit_name=unit, ...)"
    units::NamedTuple
    "sub-hydrological routing element: (route_name=route, ...)"
    routes::NamedTuple
    "The area corresponding to the node"
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
    timeidx::Vector=collect(1:length(input[first(keys(input))])),
    solver::AbstractSolver=ODESolver()
)
    output_list = []
    for idx in 1:length(node.subnames) # @threads 
        subname = node.subnames[idx]
        tmp_input = input[subname]
        tmp_unit_ouput = node.units[idx](tmp_input, pas[subname], timeidx=timeidx, solver=solver)
        tmp_route_output = node.routes[idx](tmp_unit_ouput, pas[subname], timeidx=timeidx)
        push!(output_list, getproperty(tmp_route_output, :flow) .* pas[subname][:weight])
    end
    (flow=sum(output_list),)
end

export HydroNode