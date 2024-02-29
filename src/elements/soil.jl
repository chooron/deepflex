"""
SoilWaterReservoir in Exp-Hydro
"""
function Soil_ExpHydro(; name::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        # Pet([:Temp, :Lday]),
        Tranparent([:Infiltration]),
        Evap([:SoilWater, :Pet], parameters=parameters[[:Smax]]),
        Baseflow([:SoilWater], parameters=parameters[[:Smax, :Qmax, :f]]),
        Surfaceflow([:SoilWater], parameters=parameters[[:Smax]]),
        Flow([:Baseflow, :Surfaceflow])
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(SoilWater=input[:Infiltration] .- input[:Evap] .- input[:Flow])
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
function Soil_M50(; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    et_ann = Lux.Chain(Lux.Dense(3, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 1, leakyrelu))
    q_ann = Lux.Chain(Lux.Dense(2, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 1, leakyrelu))

    funcs = [
        Tranparent([:Infiltration]),
        # ET ANN
        NNFlux([:SnowWater, :SoilWater, :Temp], [:Evap], model=et_ann, seed=42),
        # Q ANN
        NNFlux([:SoilWater, :Prcp], [:Flow], model=q_ann, seed=42),
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(SoilWater=input[:Infiltration] .-
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
function Soil_M100(; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    ann = Lux.Chain(Lux.Dense(3, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 1, leakyrelu))

    funcs = [
        Tranparent([:Melt, :Rainfall, :Temp, :SnowWater, :SoilWater, :Evap, :Lday]),
        NNFlux([:SnowWater, :SoilWater, :Temp], [:Evap], model=et_ann, seed=42)
        #TODO 还没写完
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
function Soil_GR4J(; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    funcs = [
        Saturation([:SoilWater, :Infiltration], parameters=parameters[[:x1]]),
        Evap([:SoilWater, :Pet], parameters=parameters[[:x1]]),
        Percolation([:SoilWater], parameters=parameters[[:x1]]),
        SimpleFlux([:Infiltration, :Percolation, :Saturation], [:SlowFlow, :FastFlow],
            ComponentVector{T}(),
            (i, p) -> [@.((i[:Infiltration] - i[:Saturation] + i[:Percolation]) * 0.9),
                @.((i[:Infiltration] - i[:Saturation] + i[:Percolation]) * 0.1)])
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
function Soil_HyMOD(; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    funcs = [
        Saturation([:SoilWater, :Infiltration], parameters=parameters[[:Smax, :b]]),
        Evap([:SoilWater, :Pet], parameters=parameters[[:Smax]]),
        Splitter([:Saturation], [:FastFlow, :SlowFlow],
            parameters=ComponentVector(FastFlow=parameters[:a], SlowFlow=(1 .- parameters[:a])))
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(SoilWater=input[:Infiltration] .- input[:Evap] .- input[:Saturation])
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
function Soil_XAJ(; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    tmp_func = (i, p) -> begin
        free_water, flux_in = i[:FreeWater], i[:FluxIn]
        Smax, ex = p[:Smax], p[:ex]
        tmp_re = @.(step_func(1 - free_water / Smax) * (1 - free_water / Smax))
        @.[1 - (step_func(1 - tmp_re) * (1 - tmp_re) + step_func(tmp_re - 1) * (tmp_re - 1))^ex * flux_in]
    end

    funcs = [
        Saturation([:TensionWater, :Infiltration], parameters=parameters[[:Aim, :Wmax, :a, :b]]),
        Evap([:TensionWater, :Pet], parameters=parameters[[:c, :LM]]),
        SimpleFlux([Dict(:FreeWater => :FreeWater), Dict(:Saturation => :FluxIn)], [:SurfaceFlow],
            parameters=parameters[[:Smax, :ex]], func=tmp_func),
        SimpleFlux([Dict(:FreeWater => :FreeWater), Dict(:Surfaceflow => :FluxIn)], [:InterFlow],
            parameters=parameters[[:Smax, :ex]], func=tmp_func),
        SimpleFlux([Dict(:FreeWater => :FreeWater), Dict(:Interflow => :FluxIn)], [:BaseFlow],
            parameters=parameters[[:Smax, :ex]], func=tmp_func),
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(
            TensionWater=@.(input[:Infiltration] - input[:Evap] - input[:Saturation]),
            FreeWater=@.(input[:Saturation] - input[:SurfaceFlow] - input[:InterFlow] - input[:BaseFlow])
        )
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
HBV
"""
function Soil_HBV(; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    funcs = [
        SimpleFlux([:SoilWater], [:Capillary], parameters=parameters[[:cflux, :fc]], func=(i, p) -> @.[p[:cflux] * (1 - i[:SoilWater] / p[:fc])]),
        Evap([:SoilWater, :Pet], parameters=parameters[[:lp, :fc]]),
        Recharge([:SoilWater, :Infiltration], parameters=parameters[[:fc, :β]]),
        SimpleFlux([:UpperZone], [:InterFlow], parameters=parameters[[:k0, :α]], func=(i, p) -> @.[p[:k0] * i[:UpperZone]^(1 + p[:α])]),
        SimpleFlux([:LowerZone], [:BaseFlow], parameters=parameters[[:k1]], func=(i, p) -> @.[p[:k1] * i[:LowerZone]]),
        Summation([:InterFlow, :BaseFlow], [:Flow])
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(
            SoilWater=@.(input[:Infiltration] + input[:Capillary] - input[:Evap] - input[:Recharge]),
            UpperZone=@.(input[:Recharge] - input[:Capillary] - input[:InterFlow] - parameters[:c]),
            LowerZone=@.(parameters[:c] - input[:BaseFlow])
        )
    end

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end