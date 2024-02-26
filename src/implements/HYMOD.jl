function HYMOD(; name::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    elements = [
        SoilWater_HYMOD(name="sr", parameters=parameters, init_states=init_states[[:SoilWater]]),
        RoutingStore_HYMOD(name="rs", parameters=parameters, init_states=init_states[[:F1, :F2, :F3, :S1]]),
    ]
    build_unit(name=name, elements=elements)
end