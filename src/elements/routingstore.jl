"""
RoutingStore in GR4J
"""
function RoutingStore_GR4J(; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Tranparent([:Q9]),
        Baseflow([:RoutingStore], parameters=parameters[[:x3, :γ]]),
        Recharge([:RoutingStore], parameters=parameters[[:x2, :x3, :ω]]),
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(SoilWater=input[:Q9] + input[:Recharge] - input[:Baseflow])
    end

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end