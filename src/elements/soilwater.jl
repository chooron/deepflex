"""
SoilWaterReservoir in Exp-Hydro
"""
function SoilWaterReservoir_ExpHydro(; id::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Rainfall([:Prcp, :Temp], parameters=parameters[[:Tmin]], weights=ComponentVector(SoilWater=1.0)),
        Pet([:Temp, :Lday], parameters=ComponentVector{T}(), weights=ComponentVector(SoilWater=0.0)),
        Evap([:SoilWater, :Pet], parameters=parameters[[:Smax]], weights=ComponentVector(SoilWater=-1.0)),
        Baseflow([:SoilWater], parameters=parameters[[:Smax, :Qmax, :f]], weights=ComponentVector(SoilWater=-1.0)),
        Surfaceflow([:SoilWater], parameters=parameters[[:Smax]], weights=ComponentVector(SoilWater=-1.0)),
        Flow([:Baseflow, :Surfaceflow], parameters=ComponentVector{T}(), weights=ComponentVector(SoilWater=0.0))
    ]
    ODEElement(
        id=id,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs
    )
end

"""
SoilWaterReservoir in M50
"""
function SoilWaterReservoir_M50(; id::String, parameters::Dict{Symbol,T}, init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Rainfall([:Prcp, :Temp], parameters=Dict(:Tmin => parameters[:Tmin])),
        LinearNN([:SnowWater, :SoilWater, :Temp], [:Evap], hidd_size=32, hidd_layer=1),
    ]
end