"""
SnowWaterReservoir in Exp-Hydro
"""
function SnowWaterReservoir_Exphydro(; id::String, parameters::Dict{Symbol,T}, init_states::Dict{Symbol,T}) where {T<:Number}
    funcs = [
        Snowfall([:Prcp, :Temp], parameters=Dict(:Tmin => parameters[:Tmin]))
        Melt([:SoilWater, :Temp], parameters=Dict(:Tmax => parameters[:Tmax], :Df => parameters[:Df]))
    ]
    multiplier = Dict(:SoilWater=>ComponentVector{T}(Snow=1.0, Melt=-1.0))
    ODEElement{T}(
        id=id,
        init_states=init_states,
        funcs=funcs,
        multiplier=multiplier
    )
end