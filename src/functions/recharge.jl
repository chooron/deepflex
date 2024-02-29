function Recharge(input_names::Union{Vector{Symbol},Vector{Dict{Symbol,Symbol}}},
    output_names::Vector{Symbol}=[:Recharge];
    parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        output_names,
        parameters=parameters,
        func=recharge_func
    )
end

function recharge_func(
    input::gen_namedtuple_type([:RoutingStore], T),
    parameters::gen_namedtuple_type([:x2, :x3, :ω], T)
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    routing_store = input[:RoutingStore]
    x2, x3, ω = parameters[:x2], parameters[:x3], parameters[:ω]
    [@.(x2 / (x3^ω) * routing_store^ω)]
end

function recharge_func(
    input::gen_namedtuple_type([:SoilWater, :Infiltration], T),
    parameters::gen_namedtuple_type([:fc, :β], T)
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    soil_water, infiltration = input[:SoilWater], input[:Infiltration]
    fc, β = parameters[:fc], parameters[:β]
    [@.((infiltration) * (soil_water / fc)^β)]
end