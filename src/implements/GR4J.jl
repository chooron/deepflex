function GR4J(; id::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    elements = [
        SoilWater_GR4J(id="sr", parameters=parameters, init_states=init_states[[:SoilWater]]),
        RoutingStore_GR4J(id="rs", parameters=parameters, init_states=init_states[[:RoutingStore]]),
    ]
    build_unit(id=id, elements=elements)
end