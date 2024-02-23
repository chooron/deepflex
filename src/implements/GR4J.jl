function GR4J(; name::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    elements = [
        SoilWater_GR4J(name="sr", parameters=parameters, init_states=init_states[[:SoilWater]]),
        LAGElement(name="lag", lag_time=parameters[:x4], lag_func=Dict(:Q9=>DeepFlex.unit_hydro1, :Q1=>DeepFlex.unit_hydro2)),
        RoutingStore_GR4J(name="rs", parameters=parameters, init_states=init_states[[:RoutingStore]]),
    ]
    build_unit(name=name, elements=elements)
end