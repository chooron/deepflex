"""
Exp-Hydro model
"""
function ExpHydro(; id::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    elements = [
        SnowWaterReservoir_ExpHydro(id="sr", parameters=parameters, init_states=init_states),
        SoilWaterReservoir_ExpHydro(id="wr", parameters=parameters, init_states=init_states)
    ]
    build_unit(id=id, elements=elements)
end
