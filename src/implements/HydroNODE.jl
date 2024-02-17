"""
Implement for [Improving hydrologic models for predictions and process understanding using neural ODEs](https://hess.copernicus.org/articles/26/5085/2022/)
"""
function M50(; id::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    elements = [
        SnowWater_ExpHydro(id="sr", parameters=parameters, init_states=init_states[[:SnowWater]]),
        SoilWater_M50(id="wr", parameters=parameters, init_states=init_states[[:SoilWater]])
    ]
    build_unit(id=id, elements=elements)
end

function M100(; id::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    elements = [
        LinearNN([:SnowWater, :SoilWater, :Temp, :Prcp], [:Snowfall, :Rainfall, :Melt, :Evap, :Flow], hidd_size=32, hidd_layer=1),
        SnowWater_ExpHydro_ODE(id="sr", parameters=parameters, init_states=init_states[[:SnowWater]]),
        SoilWater_M50_ODE(id="wr", parameters=parameters, init_states=init_states[[:SoilWater]])
    ]
    build_unit(id=id, elements=elements)
end