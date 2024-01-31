"""
Implement for [Improving hydrologic models for predictions and process understanding using neural ODEs](https://hess.copernicus.org/articles/26/5085/2022/)
"""
function M50(; id::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    elements = [
        SnowWater_ExpHydro_ODE(id="sr", parameters=parameters, init_states=init_states[[:SnowWater]]),
        SoilWater_M50_ODE(id="wr", parameters=parameters, init_states=init_states[[:SoilWater]])
    ]
    build_unit(id=id, elements=elements)
end