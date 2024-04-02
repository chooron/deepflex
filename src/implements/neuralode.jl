"""
Implement for [Improving hydrologic models for predictions and process understanding using neural ODEs](https://hess.copernicus.org/articles/26/5085/2022/)
"""
function M50(; name::Symbol)
    elements = [
        Surface_ExpHydro(name=name),
        Soil_M50(name=name)
    ]
    HydroUnit(name, elements=elements)
end

function M100(; name::Symbol, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    elements = [
        LinearNN(
            [:SnowWater, :SoilWater, :Temp, :Prcp],
            [:Snowfall, :Rainfall, :Melt, :Evap, :Flow],
            hnamed_size=32,
            hnamed_layer=1
        ),
        SnowWater_ExpHydro_ODE(
            name=:sr,
            parameters=parameters,
            init_states=init_states[[:SnowWater]]
        ),
        SoilWater_M50_ODE(
            name=:wr,
            parameters=parameters,
            init_states=init_states[[:SoilWater]]
        )
    ]
    build_unit(name=name, elements=elements)
end