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
    i::gen_namedtuple_type([:baseflow, :surfaceflow], T),
    p::NamedTuple,
    sf::Function
)::Union{T,Vector{T}} where {T<:Number}
    i[:baseflow] .+ i[:surfaceflow]
end

function flow_func(
    i::gen_namedtuple_type([:routedflow, :recharge, :fastflow], T),
    p::NamedTuple,
    sf::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(i[:routedflow] + sf(i[:fastflow] + i[:recharge]) * (i[:fastflow] + i[:recharge]))
end

function flow_func(
    i::gen_namedtuple_type([:baseflow, :interflow, :surfaceflow], T),
    p::NamedTuple,
    sf::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(i[:surfaceflow] + i[:baseflow] + i[:interflow])
end
