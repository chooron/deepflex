"""
SnowWaterReservoir in Exp-Hydro
"""
function Surface_ExpHydro(; name::Symbol,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    funcs = [
        Pet([:Temp, :Lday]),
        Snowfall([:Prcp, :Temp], parameters=parameters[[:Tmin]]),
        Melt([:SnowWater, :Temp], parameters=parameters[[:Tmax, :Df]]),
        Rainfall([:Prcp, :Temp], parameters=parameters[[:Tmin]]),
        Infiltration([:Rainfall, :Melt])
    ]

    d_funcs = [
        Differ(Dict(:In => [:Snowfall], :Out => [:Melt]), [:SnowWater]),
    ]

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
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
        SimpleFlux([:Prcp, :Pet], [:Pet], parameters=ComponentVector(),
            func=(i, p, sf) -> @.(sf(i[:Pet] - i[:Prcp]) * (i[:Pet] - i[:Prcp]))),
        Infiltration([:Rainfall])
    ]

    SimpleElement(
        name=name,
        parameters=ComponentVector(),
        funcs=funcs
    )
end


function Surface_HBV(; name::Symbol,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    funcs = [
        Snowfall([:Prcp, :Temp], parameters=parameters[[:tt, :tti]]),
        SimpleFlux([:Temp], [:Refreeze], parameters=parameters[[:cfr, :cfmax, :ttm]],
            func=(i, p, sf) -> @.(sf(p[:ttm] - i[:Temp]) * p[:cfr] * p[:cfmax] * (p[:ttm] - i[:Temp]))),
        Melt([:Temp], parameters=parameters[[:cfmax, :ttm]]),
        Rainfall([:Prcp, :Temp], parameters=parameters[[:tt, :tti]]),
        Infiltration([:SnowWater, :LiquidWater, :Rainfall, :Melt], parameters=parameters[[:whc]]),
    ]

    d_funcs = [
        Differ(Dict(:In => [:Snowfall, :Refreeze], :Out => [:Melt]), [:SnowWater]),
        Differ(Dict(:In => [:Rainfall, :Melt], :Out => [:Refreeze, :Infiltration]), [:LiquidWater]),
    ]

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        d_funcs=d_funcs
    )
end

function Surface_XAJ(; name::Symbol,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    funcs = [
        Snowfall([:Prcp, :Temp], parameters=parameters[[:tt, :tti]]),
        SimpleFlux([:Temp], [:Refreeze], parameters=parameters[[:cfr, :cfmax, :ttm]],
            func=(i, p, sf) -> @.(sf(p[:ttm] - i[:Temp]) * p[:cfr] * p[:cfmax] * (p[:ttm] - i[:Temp]))),
        Melt([:Temp], parameters=parameters[[:cfmax, :ttm]]),
        Rainfall([:Prcp, :Temp], parameters=parameters[[:tt, :tti]]),
        Infiltration([:LiquidWater, :Rainfall, :Melt], parameters=parameters[[:whc, :sp]]),
    ]

    d_funcs = [
        Differ(Dict(:In => [:Snowfall, :Refreeze], :Out => [:Melt]), [:SnowWater]),
        Differ(Dict(:In => [:Rainfall, :Melt], :Out => [:Refreeze, :Infiltration]), [:LiquidWater]),
    ]

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end


function Surface_M100(; name::Symbol,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    funcs = [
        Tranparent([:Melt, :Snowfall, :Temp, :SnowWater]),
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(SnowWater=begin
            relu(sinh(input[:Snowfall]) .* step_fct(input[:Temp]))
            .-relu(step_func(input[:SnowWater]) .* sinh(input[:Melt]))
        end)
    end

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end