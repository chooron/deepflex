function GR4J(; name::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    elements = [
        Surface_GR4J(name="sf"),
        Soil_GR4J(name="sl", parameters=parameters[[:x1]], init_states=init_states[[:SoilWater]]),
        Lag_GR4J(name="lg", parameters=parameters[[:x4]]),
        Routing_GR4J(name="rt", parameters=parameters[[:x2, :x3, :ω, :γ]], init_states=init_states[[:RoutingStore]]),
    ]
    build_unit(name=name, elements=elements)
end