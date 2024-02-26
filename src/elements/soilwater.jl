"""
SoilWaterReservoir in Exp-Hydro
"""
function SoilWater_ExpHydro(; name::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Rainfall([:Prcp, :Temp], parameters=parameters[[:Tmin]]),
        Tranparent([:Melt]),
        Pet([:Temp, :Lday]),
        Evap([:SoilWater, :Pet], parameters=parameters[[:Smax]]),
        Baseflow([:SoilWater], parameters=parameters[[:Smax, :Qmax, :f]]),
        Surfaceflow([:SoilWater], parameters=parameters[[:Smax]]),
        Flow([:Baseflow, :Surfaceflow])
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(SoilWater=input[:Rainfall] .+ input[:Melt] .- input[:Evap] .- input[:Flow])
    end

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end


"""
SoilWaterReservoir in M50
"""
function SoilWater_M50(; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    et_ann = Lux.Chain(Lux.Dense(3, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 1, leakyrelu))
    q_ann = Lux.Chain(Lux.Dense(2, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 1, leakyrelu))

    funcs = [
        Rainfall([:Prcp, :Temp], parameters=parameters[[:Tmin]]),
        # ET ANN
        NNFlux([:SnowWater, :SoilWater, :Temp], [:Evap], model=et_ann, seed=42),
        # Q ANN
        NNFlux([:SoilWater, :Prcp], [:Flow], model=q_ann, seed=42),
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(SoilWater=input[:Rainfall] .+
                                  input[:Melt] .-
                                  step_func(input[:SoilWater]) .* input[:Lday] .* exp(input[:Evap]) .-
                                  step_func(input[:SoilWater]) .* exp(input[:Flow]))
    end

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end

"""
SoilWaterReservoir in M100
"""
function SoilWater_M100(; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Tranparent([:Melt, :Rainfall, :Temp, :SnowWater, :SoilWater, :Evap, :Lday]),
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(SoilWater=begin
            relu(sinh(input[:Rainfall])) .+
            relu(step_func(input[:SnowWater]) * sinh(input[:Melt])) .-
            step_fct(input[:SoilWater]) .* input[:Lday] .* exp(input[:Evap]) .-
            step_fct(input[:SoilWater]) .* exp(input[:Flow])
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

"""
SoilWaterReservoir in GR4J
"""
function SoilWater_GR4J(; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    funcs = [
        Rainfall([:Prcp, :Pet]),
        Saturation([:SoilWater, :Rainfall], parameters=parameters[[:x1]]),
        Evap([:SoilWater, :Prcp, :Pet], parameters=parameters[[:x1]]),
        Percolation([:SoilWater], parameters=parameters[[:x1]]),
        SimpleFlux([:Rainfall, :Percolation, :Saturation], [:Q9, :Q1],
            ComponentVector{T}(),
            (i, p) -> [@.((i[:Rainfall] - i[:Saturation] + i[:Percolation]) * 0.9),
                @.((i[:Rainfall] - i[:Saturation] + i[:Percolation]) * 0.1)])
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(SoilWater=input[:Rainfall] .- input[:Evap] .- input[:Percolation])
    end

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end

"""
SoilWaterReservoir in HYMOD
"""
function SoilWater_HYMOD(; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    funcs = [
        Rainfall([:Prcp]),
        Saturation([:SoilWater, :Rainfall], parameters=parameters[[:Smax, :b]]),
        Evap([:SoilWater, :Pet], parameters=parameters[[:Smax]])
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(SoilWater=input[:Rainfall] .- input[:Evap] .- input[:Saturation])
    end

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end


"""
SoilWaterReservoir in XAJ
"""
function TensionWater_XAJ(; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    funcs = [
        Saturation([:TensionWater, :Prcp], parameters=parameters[[:Aim, :Wmax, :a, :b]]),
        Evap([:TensionWater, :Pet], parameters=parameters[[:c, :LM]])
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(TensionWater=@.(input[:Prcp] * (1 - parameters[:Aim]) - input[:Evap] - input[:Saturation]))
    end

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end

function FreeWater_XAJ(; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    funcs = [
        Saturation([Dict(:FreeWater => :FreeWater), Dict(:Saturation => :FluxIn)], [:SurfaceRunoff], parameters=parameters[[:Smax, :ex]]),
        Saturation([Dict(:FreeWater => :FreeWater), Dict(:SurfaceRunoff => :FluxIn)], [:InterRunoff], parameters=parameters[[:Smax, :ex]]),
        Saturation([Dict(:FreeWater => :FreeWater), Dict(:InterRunoff => :FluxIn)], [:BaseRunoff], parameters=parameters[[:Smax, :ex]]),
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(FreeWater=@.(input[:Saturation] - input[:SurfaceRunoff] - input[:InnerRunoff] - input[:BaseRunoff]))
    end

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end