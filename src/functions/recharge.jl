function Recharge(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux{T}(
        input_names,
        [:Recharge],
        parameters,
        recharge_func
    )
end

function recharge_func(
    input::ComponentVector{T,Vector{T},Tuple{Axis{(RoutingStore=1,)}}},
    parameters::ComponentVector{T,Vector{T},Tuple{Axis{(x2=1, x3=2, ω=3)}}}
) where {T<:Number}
    routing_store = input[:RoutingStore]
    x2, x3, ω = parameters[:x2], parameters[:x3], parameters[:ω]
    ComponentVector(Recharge=x2 / (x3^ω) * routing_store^ω)
end