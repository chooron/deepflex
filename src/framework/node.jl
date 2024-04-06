# info HydroNode一些属性在单独计算时是用不上的，比如面积
struct HydroNode <: AbstractComponent
    name::Symbol
    units::Vector{HydroUnit}
    weight::NamedTuple
    route::NamedTuple
end

# function Node(name::Symbol; units::Vector{HydroUnit}, weight::NamedTuple, route::NamedTuple)
#     Node(name, units, weight, route)
# end

function (node::Node)(input::NamedTuple, params::NamedTuple, init_states::NamedTuple)
    unit_names = [u.name for u in node.units]
    unit_ouputs = namedtuple(unit_names, [
        begin
            tmp_ouput = u(input, params, init_states)
            route_output = node.route[u.name](tmp_ouput, params, init_states) .* node.weight[u.name]
            route_output
        end
        for u in node.units
    ])
    unit_ouputs
end