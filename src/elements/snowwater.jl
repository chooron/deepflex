"""
SnowWaterReservoir in Exp-Hydro
"""
function SnowWaterReservoir_ExpHydro(; id::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Snowfall([:Prcp, :Temp], parameters=parameters[[:Tmin]], weights=ComponentVector(Snowfall=1.0))
        Melt([:SnowWater, :Temp], parameters=parameters[[:Tmax, :Df]], weights=ComponentVector(Melt=-1.0))
    ]
    ODEElement(
        id=id,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs
    )
end