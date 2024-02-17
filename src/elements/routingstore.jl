"""
RoutingStore in GR4J
"""
function RoutingStore_GR4J(; id::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Splitter([:Flow], parameters=ComponentVector{T}(Q9=0.9, Q1=0.1)),
        Routing([:Q9], lag_time=parameters[:x4], lag_func=unit_hydro1),
        Routing([:Q1], lag_time=parameters[:x4], lag_func=unit_hydro2),
        Baseflow([:RoutingStore], parameters=parameters[[:x3, :γ]]),
        Recharge([:RoutingStore], parameters=parameters[[:x2, :x3, :ω]]),
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(SoilWater=input[:Q9] + input[:Recharge] - input[:Baseflow])
    end

    build_element(
        id=id,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end