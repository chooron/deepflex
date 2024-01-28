"""
SoilWaterReservoir in Exp-Hydro
"""
function SoilWaterReservoir_ExpHydro(; id::String, parameters::Dict{Symbol,T}, init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Rainfall([:Prcp, :Temp], parameters=Dict(:Tmin => parameters[:Tmin])),
        Pet([:Temp, :Lday], parameters=Dict()),
        Evap([:SoilWater, :Pet], parameters=Dict(:Smax => parameters[:Smax])),
        Baseflow([:SoilWater], parameters=Dict(:Smax => parameters[:Smax], :Qmax => parameters[:Qmax], :f => parameters[:f])),
        Surfaceflow([:SoilWater], parameters=Dict(:Smax => parameters[:Smax])),
        Flow([:Baseflow, :Surfaceflow], parameters=Dict())
    ]
    multiplier = Dict(:SoilWater => ComponentVector{T}(Rain=1.0, Melt=1.0, Evap=-1.0, Qb=-1.0, Qs=-1.0))
    ODEElement{T}(
        id=id,
        init_states=init_states,
        funcs=funcs,
        multiplier=multiplier
    )
end

"""
SoilWaterReservoir in M50
"""
function SoilWaterReservoir_M50(; id::String, parameters::Dict{Symbol,T}, init_states::ComponentVector{T})
    funcs = [
        Rainfall([:Prcp, :Temp], parameters=Dict(:Tmin => parameters[:Tmin])),
        LinearNN([:SnowWater, :SoilWater, :Temp], [:Evap], hidd_size=32, hidd_layer=1),
    ]
end