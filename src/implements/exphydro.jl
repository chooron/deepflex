"""
Exp-Hydro model
"""
function ExpHydro(; name::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    elements = [
        SnowWater_ExpHydro(name="sr", parameters=parameters, init_states=init_states[[:SnowWater]]),
        SoilWater_ExpHydro(name="wr", parameters=parameters, init_states=init_states[[:SoilWater]])
    ]
    build_unit(name=name, elements=elements)
end
