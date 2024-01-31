"""
SoilWaterReservoir in Exp-Hydro
"""
function SoilWater_ExpHydro_ODE(; id::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Rainfall([:Prcp, :Temp], parameters=parameters[[:Tmin]]),
        Tranparent([:Melt], parameters=ComponentVector{T}()),
        Pet([:Temp, :Lday], parameters=ComponentVector{T}()),
        Evap([:SoilWater, :Pet], parameters=parameters[[:Smax]]),
        Baseflow([:SoilWater], parameters=parameters[[:Smax, :Qmax, :f]]),
        Surfaceflow([:SoilWater], parameters=parameters[[:Smax]]),
        Flow([:Baseflow, :Surfaceflow], parameters=ComponentVector{T}())
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(SoilWater=input[:Rainfall] + input[:Melt] - input[:Evap] - input[:Flow])
    end

    ODEElement(
        id=id,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end

"""
SoilWaterReservoir in M50
"""
function SoilWater_M50_ODE(; id::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Rainfall([:Prcp, :Temp], parameters=parameters[[:Tmin]]),
        # ET ANN
        LinearNN([:SnowWater, :SoilWater, :Temp], [:Evap], hidd_size=32, hidd_layer=1),
        # Q ANN
        LinearNN([:SoilWater, :Prcp], [:Flow], hidd_size=32, hidd_layer=1),
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(SoilWater=input[:Rainfall] + input[:Melt] - step_func(input[:SoilWater]) * input[:Lday] * exp(input[:Evap]) - step_func(input[:SoilWater]) * exp(input[:Flow]))
    end

    ODEElement(
        id=id,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end