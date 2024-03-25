function FlowFlux(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:flow;
    param_names::Vector{Symbol}=Symbol[])
    SimpleFlux(
        input_names,
        output_names,
        param_names=param_names,
        func=flow_func,
    )
end

function flow_func(
    i::namedtuple(:baseflow, :surfaceflow),
    p::NamedTuple,
    sf::Function
)
    i[:baseflow] .+ i[:surfaceflow]
end

function flow_func(
    i::namedtuple(:routedflow, :recharge, :fastflow),
    p::NamedTuple,
    sf::Function
)
    @.(i[:routedflow] + sf(i[:fastflow] + i[:recharge]) * (i[:fastflow] + i[:recharge]))
end

function flow_func(
    i::namedtuple(:baseflow, :interflow, :surfaceflow),
    p::NamedTuple,
    sf::Function
)
    @.(i[:surfaceflow] + i[:baseflow] + i[:interflow])
end
