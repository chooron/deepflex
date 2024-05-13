function FlowFlux(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:flow;
    param_names::Vector{Symbol}=Symbol[],
    smooth_func::Function=step_func)
    SimpleFlux(
        input_names,
        output_names,
        param_names=param_names,
        func=flow_func,
        smooth_func=smooth_func
    )
end

function flow_func(
    i::namedtuple(:baseflow, :surfaceflow),
    p::NamedTuple;
    kw...
)
    i[:baseflow] .+ i[:surfaceflow]
end

function flow_func(
    i::namedtuple(:routedflow, :recharge, :fastflow),
    p::NamedTuple;
    kw...
)
    @.(i[:routedflow] + sf(i[:fastflow] + i[:recharge]) * (i[:fastflow] + i[:recharge]))
end

function flow_func(
    i::namedtuple(:baseflow, :interflow, :surfaceflow),
    p::NamedTuple;
    kw...
)
    @.(i[:surfaceflow] + i[:baseflow] + i[:interflow])
end

export FlowFlux, flow_func