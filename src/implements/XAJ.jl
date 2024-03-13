function XAJ(; name::Symbol, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    elements = [
        Surface_GR4J(name=:sf),
        Soil_XAJ(name=:sl,
            parameters=parameters[[:Aim, :Wmax, :a, :b, :c, :LM, :Smax, :ex]],
            init_states=init_states[[:TensionWater, :FreeWater]]),
        Routing_XAJ(name=:rt,
            parameters=parameters[[:ci, :cg, :Aim]],
            init_states=init_states[[:InterRouting, :BaseRouting]])
    ]
    build_unit(name=name, elements=elements)
end