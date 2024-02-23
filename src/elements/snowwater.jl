"""
SnowWaterReservoir in Exp-Hydro
"""
function SnowWater_ExpHydro(; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    funcs = [
        Snowfall([:Prcp, :Temp], parameters=parameters[[:Tmin]])
        Melt([:SnowWater, :Temp], parameters=parameters[[:Tmax, :Df]])
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(SnowWater=input[:Snowfall] .- input[:Melt],)
    end

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end

function SnowWater_M100(; name::String,
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

function SnowWater_M100(; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}
    
    funcs = [
        Tranparent([:Melt, :Snowfall, :Temp, :SnowWater]),
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(SnowWater=begin
            relu(sinh(input[:Snowfall]) * step_fct(input[:Temp]))
            -relu(step_func(input[:SnowWater]) * sinh(input[:Melt]))
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