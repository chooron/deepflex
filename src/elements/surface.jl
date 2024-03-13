"""
SnowWaterReservoir in Exp-Hydro
"""
function Surface_ExpHydro(; name::Symbol)
    funcs = [
        Pet([:Temp, :Lday]),
        Snowfall([:Prcp, :Temp], parameter_names=[:Tmin]),
        Melt([:SnowWater, :Temp], parameter_names=[:Tmax, :Df]),
        Rainfall([:Prcp, :Temp], parameter_names=[:Tmin]),
        Infiltration([:Rainfall, :Melt])
    ]

    d_funcs = [
        Differ(Dict(:In => [:Snowfall], :Out => [:Melt]), [:SnowWater]),
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        d_funcs=d_funcs
    )
end

"""
SnowWaterReservoir in Exp-Hydro
"""
function Surface_GR4J(; name::Symbol)
    funcs = [
        Rainfall([:Prcp, :Pet]),
        SimpleFlux([:Prcp, :Pet], :Pet, parameters=Symbol[],
            func=(i, p, sf) -> @.(sf(i[:Pet] - i[:Prcp]) * (i[:Pet] - i[:Prcp]))),
        Infiltration([:Rainfall])
    ]

    SimpleElement(
        name=name,
        funcs=funcs
    )
end


function Surface_HBV(; name::Symbol)
    funcs = [
        Snowfall([:Prcp, :Temp], parameter_names=[:tt, :tti]),
        SimpleFlux([:Temp], :Refreeze, parameter_names=[:cfr, :cfmax, :ttm],
            func=(i, p, sf) -> @.(sf(p[:ttm] - i[:Temp]) * p[:cfr] * p[:cfmax] * (p[:ttm] - i[:Temp]))),
        Melt([:Temp], parameter_names=[:cfmax, :ttm]),
        Rainfall([:Prcp, :Temp], parameter_names=[:tt, :tti]),
        Infiltration([:SnowWater, :LiquidWater, :Rainfall, :Melt], parameter_names=[:whc]),
    ]

    d_funcs = [
        Differ(Dict(:In => [:Snowfall, :Refreeze], :Out => [:Melt]), [:SnowWater]),
        Differ(Dict(:In => [:Rainfall, :Melt], :Out => [:Refreeze, :Infiltration]), [:LiquidWater]),
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        d_funcs=d_funcs
    )
end

function Surface_XAJ(; name::Symbol)
    funcs = [
        Snowfall([:Prcp, :Temp], parameter_names=[:tt, :tti]),
        SimpleFlux([:Temp], :Refreeze, parameter_names=[:cfr, :cfmax, :ttm],
            func=(i, p, sf) -> @.(sf(p[:ttm] - i[:Temp]) * p[:cfr] * p[:cfmax] * (p[:ttm] - i[:Temp]))),
        Melt([:Temp], parameter_names=[:cfmax, :ttm]),
        Rainfall([:Prcp, :Temp], parameter_names=[:tt, :tti]),
        Infiltration([:LiquidWater, :Rainfall, :Melt], parameter_names=[:whc, :sp]),
    ]

    d_funcs = [
        Differ(Dict(:In => [:Snowfall, :Refreeze], :Out => [:Melt]), [:SnowWater]),
        Differ(Dict(:In => [:Rainfall, :Melt], :Out => [:Refreeze, :Infiltration]), [:LiquidWater]),
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        d_funcs=d_funcs
    )
end


function Surface_M100(; name::Symbol)
    funcs = SimpleFlux[]

    d_funcs = [
        SimpleFlux([:SnowWater, :Snowfall, :Temp, :Melt], :SnowWater, parameters=Symbol[],
            func=(i, p, sf) -> @.(relu(sinh(i[:Snowfall]) * sf(i[:Temp])) - relu(sf(i[:SnowWater]) * sinh(i[:Melt])))),
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        d_funcs=d_funcs
    )
end