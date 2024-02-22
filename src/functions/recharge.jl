function Recharge(input_names::Vector{Symbol}; parameters::Union{ComponentVector{T},Nothing}=nothing) where {T<:Number}
    build_flux(
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

function recharge_func(
    input::gen_namedtuple_type([:RoutingStore], T),
    parameters::gen_namedtuple_type([:x2,:x3,:ω], T)
)::gen_namedtuple_type([:Recharge], T) where {T<:Number}
    routing_store = input[:RoutingStore]
    x2, x3, ω = parameters[:x2], parameters[:x3], parameters[:ω]
    (Recharge=@.(x2 / (x3^ω) * routing_store^ω),)
end