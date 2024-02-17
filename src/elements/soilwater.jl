"""
SoilWaterReservoir in Exp-Hydro
"""
function SoilWater_ExpHydro(; id::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
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

    build_element(
        id=id,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end

function SoilWater_ExpHydro_MTK(; id::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}


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

    build_element(
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
function SoilWater_M50(; id::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
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

    build_element(
        id=id,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end

"""
SoilWaterReservoir in M100
"""
function SoilWater_M100(; id::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Tranparent([:Melt, :Rainfall, :Temp, :SnowWater, :SoilWater, :Evap, :Lday], parameters=ComponentVector{T}()),
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(SoilWater=begin
            relu(sinh(input[:Rainfall])) +
            relu(step_func(input[:SnowWater]) * sinh(input[:Melt]))
            -step_fct(input[:SoilWater]) * input[:Lday] * exp(input[:Evap])
            -step_fct(input[:SoilWater]) * exp(input[:Flow])
        end)
    end

    build_element(
        id=id,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end

"""
SoilWaterReservoir in GR4J
"""
function SoilWater_GR4J(; id::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Rainfall([:Prcp, :Pet], parameters=parameters[[:x1]])
        Saturation([:SoilWater, :Rainfall, :Pet], parameters=parameters[[:x1]])
        Evap([:SoilWater, :Prcp, :Pet], parameters=parameters[[:x1]])
        Percolation([:SoilWater], parameters=parameters[[:x1]])
        Flow([:Rainfall, :Saturation, :7], parameters=ComponentVector{T}())
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(SoilWater=input[:Rainfall] - input[:Evap] - input[:Percolation])
    end

    build_element(
        id=id,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end