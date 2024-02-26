function HYMOD(; name::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    elements = [
        TensionWater_XAJ(name="tw", parameters=parameters, init_states=init_states[[:TensionWater]]),
        FreeWater_XAJ(name="fw", parameters=parameters, init_states=init_states[[:FreeWater]]),
        RoutingStore_XAJ(name="rs", parameters=parameters, init_states=init_states[[:InterRoutingStore, :BaseRoutingStore]]),
        SimpleElement(name="out", parameters=parameters, funcs=[SimpleFlux([:Surfaceflow, :Interflow, :Baseflow], [:Flow], parameters[:Aim])])
    ]
    build_unit(name=name, elements=elements)
end