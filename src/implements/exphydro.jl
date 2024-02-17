"""
Exp-Hydro model
"""
function ExpHydro(; id::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    elements = [
        SnowWater_ExpHydro(id="sr", parameters=parameters, init_states=init_states[[:SnowWater]]),
        SoilWater_ExpHydro(id="wr", parameters=parameters, init_states=init_states[[:SoilWater]])
    ]
    build_unit(id=id, elements=elements)
end
