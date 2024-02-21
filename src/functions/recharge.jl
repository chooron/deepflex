function Recharge(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    SimpleFlux(
        input_names,
        [:Recharge],
        parameters,
        recharge_func
    )
end

function recharge_func(
    input::(@NamedTuple{RoutingStore::Union{T,Vector{T}}}),
    parameters::(@NamedTuple{x2::Union{T,Vector{T}}, x3::Union{T,Vector{T}}, ω::Union{T,Vector{T}}})
)::(@NamedTuple{Recharge::Union{T,Vector{T}}}) where {T<:Number}
    routing_store = input[:RoutingStore]
    x2, x3, ω = parameters[:x2], parameters[:x3], parameters[:ω]
    (Recharge=@.(x2 / (x3^ω) * routing_store^ω),)
end