"""
SnowWaterReservoir in Exp-Hydro
"""
function SnowWaterReservoir_ExpHydro(; id::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Snowfall([:Prcp, :Temp], parameters=parameters[[:Tmin]], weights=ComponentVector(SnowWater=1.0))
        Melt([:SnowWater, :Temp], parameters=parameters[[:Tmax, :Df]], weights=ComponentVector(SnowWater=-1.0))
    ]
    ODEElement(
        id=id,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs
    )
end