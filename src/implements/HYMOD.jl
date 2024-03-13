function HyMOD(; name::Symbol, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    elements = [
        Surface_GR4J(name="sf"),
        Soil_HyMOD(name="sl", parameters=parameters[[:Smax, :a, :b]], init_states=init_states[[:SoilWater]]),
        Routing_HyMOD(name="rt", parameters=parameters[[:kf, :ks]],
            init_states=init_states[[:FastRouting1, :FastRouting2, :FastRouting3, :SlowRouting]]),
    ]
    build_unit(name=name, elements=elements)
end