"""
SnowWaterReservoir in Exp-Hydro
"""
function SnowWater_ExpHydro_ODE(; id::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Snowfall([:Prcp, :Temp], parameters=parameters[[:Tmin]])
        Melt([:SnowWater, :Temp], parameters=parameters[[:Tmax, :Df]])
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(SnowWater=input[:Snowfall] - input[:Melt])
    end

    ODEElement(
        id=id,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end